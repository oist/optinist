from fastapi import APIRouter, WebSocket
from pydantic import BaseModel
from typing import Dict, List, Optional
import traceback
import os
import sys
import json
import yaml
import subprocess

from wrappers import wrapper_dict
from collections import OrderedDict
from wrappers.data_wrapper import *
from wrappers.optinist_exception import AlgorithmException
from .utils import get_algo_network, run_algorithm, get_results

import time
sys.path.append('../../optinist')

router = APIRouter()



def current_time():
    print(time.perf_counter())

def create_snakefile_config_from_flowlist(flowList):
    '''
    flowListを受け取り、Snakemakeに渡すconfig.yamlとして出力する。
    '''
    nodeList = flowList.get('nodeList')

    flow_config = {}
    rules_to_execute = {}
    prev_algo_output = None
    for i, item in enumerate(nodeList):
        if item["data"]["type"] == 'image':
            initial_input = "/app/" + item["data"]["path"]
        elif item["data"]["type"] == 'algo':
            algo_name = item["data"]["label"]
            if i == 1:
                algo_input = initial_input
            else:
                algo_input = prev_algo_output
            
            output_base_path = f"/app/files/{algo_name}"
            
            if not os.path.exists(output_base_path):
                print(f"Creating {output_base_path}")
                os.makedirs(output_base_path)
            
            algo_output = os.path.join(output_base_path, f"{algo_name}_out.pkl")
            
            rules_to_execute[algo_name] = {   
                    "rule_file": f"rules/{algo_name}.smk",
                    "input": algo_input,
                    "param": item["data"]["param"],
                    "output": algo_output
                }

            # 次のalgoのinputに渡すため、現在の出力ファイルのパスを保存
            prev_algo_output = algo_output

    flow_config["rules"] = rules_to_execute
    flow_config["last_output"] = algo_output

    with open('./workflow/config/config.yaml', "w") as f:
        yaml.dump(flow_config, f)

def run_snakemake(flowList, use_conda=False):
    # flowListの内容をconfig.yamlに出力
    create_snakefile_config_from_flowlist(flowList)

    # snakemake実行
    print("Executing Snakefile")
    if use_conda:
        result = subprocess.run("snakemake --cores 1 --use-conda", shell=True, cwd="./workflow")
    else:
        result = subprocess.run("snakemake --cores 1", shell=True, cwd="./workflow")

    # 以下report.html生成用（対応要検討）
    # subprocess.run("snakemake --cores 1 --report report.html", shell=True, cwd="./workflow")
    return result

def get_flow_outputs(flowList, flow_config: dict):
    flow_rules = flow_config.get("rules")
    
    flowoutputs = dict()
    
    for flowitem in flowList:
        if flowitem.label == "ImageData":
            continue
        else:
            output_pkl = flow_rules.get(flowitem.label).get("output")
            flowoutputs[flowitem.label] = output_pkl

    return flowoutputs


@router.get("/run/ready/{request_id}")
async  def run_ready(request_id: str):
    return {'message': 'ready...', 'status': 'ready', "requestId": request_id}


@router.websocket("/run")
async def websocket_endpoint(websocket: WebSocket):
    # import run_pipeline

    await websocket.accept()
    # Wait for any message from the client
    flowList = await websocket.receive_text()
    # flowList = list(map(lambda x: FlowItem(**x), json.loads(flowList)))
    flowList = json.loads(flowList)

    # snakemake実行
    run_snakemake(flowList, use_conda=False)

    graph, startNodeList, nodeDict = get_algo_network(flowList)

    try:
        item = nodeDict[startNodeList[0]]
        info = None
        prev_info = None

        while True:

            await websocket.send_json({
                'message': item['data']['label'] + ' started',
                'status': 'ready'
            })

            # run algorithm
            info = run_algorithm(info, prev_info, item)
            prev_info = info

            assert info is not None

            results = get_results(info, item)

            # Send message to the client
            await websocket.send_json({
                'message': item['data']['label'] + ' success',
                'status': 'success',
                'outputPaths': results
            })

            if len(graph[item['id']]) > 0:
                node_id = graph[item['id']][0]
                item = nodeDict[node_id]
            else:
                break

        await websocket.send_json({'message': 'completed', 'status': 'completed'})

    except AlgorithmException as e:
        await websocket.send_json({'message': e.get_message(), 'name': item['data']['label'], 'status': 'error'})
    except Exception as e:
        message  = list(traceback.TracebackException.from_exception(e).format())[-1]
        print(traceback.format_exc())
        await websocket.send_json({'message': message, 'name': item['data']['label'], 'status': 'error'})
    finally:
        print('Bye..')
        await websocket.close()
