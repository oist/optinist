from fastapi import APIRouter, WebSocket
from pydantic import BaseModel
from typing import Dict, List, Optional
import traceback
import sys
import json
from wrappers import wrapper_dict
from collections import OrderedDict
from wrappers.data_wrapper import *
from wrappers.optinist_exception import AlgorithmException

sys.path.append('../../optinist')

router = APIRouter()


class FlowItem(BaseModel):
    label: str
    path: Optional[str] = None
    type: str
    param: Optional[Dict] = {}


def get_dict_leaf_value(root_dict: dict, path_list: List[str]):
    path = path_list.pop(0)
    if(len(path_list) > 0):
        return get_dict_leaf_value(root_dict[path],path_list)
    else:
        return root_dict[path]


@router.get("/run/ready/{request_id}")
async  def run_ready(request_id: str):
    return {'message': 'ready...', 'status': 'ready', "requestId": request_id}


@router.websocket("/run")
async def websocket_endpoint(websocket: WebSocket):
    # import run_pipeline

    await websocket.accept()
    # Wait for any message from the client
    flowList = await websocket.receive_text()
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


def get_algo_network(flowList):
    flowList = json.loads(flowList)
    nodeDict = {}
    startNodeList = []

    edgeList = flowList['edgeList']
    graph = {}

    # nodeを初期化
    for node in flowList['nodeList']:
        nodeDict[node['id']] = node
        graph[node['id']] = []
        if node['type'] != 'AlgorithmNode':
            startNodeList.append(node['id'])

    # 隣接リストを登録
    for edge in edgeList:
        graph[edge['source']].append(edge['target'])

    return graph, startNodeList, nodeDict


def run_algorithm(info, prev_info, item):
    if item['type'] == 'ImageFileNode':
        info = {'path': ImageData(item['data']['path'], '')}
    elif item['type'] == 'AlgorithmNode':
        # parameterをint, floatに変換
        from .utils import string_to_float
        params = string_to_float(item['data']['param'])

        wrapper = get_dict_leaf_value(
            wrapper_dict, item['data']['path'].split('/'))
        info = wrapper["function"](
            *prev_info.values(), params=params)
    else:
        assert False, 'run_algorithm error'

    return info


def get_results(info, item):
    results = OrderedDict()
    path = item['data']['path']
    results[path] = {}

    for k, v in info.items():
        if type(v) is ImageData:
            print("ImageData")
            results[path][k] = {}
            results[path][k]['path'] = v.json_path
            results[path][k]['type'] = 'images'
            results[path][k]['max_index'] = len(v.data)
        elif type(v) is TimeSeriesData:
            print("TimeSeriesData")
            results[path][k] = {}
            results[path][k]['path'] = v.path
            results[path][k]['type'] = 'timeseries'
        elif type(v) is CorrelationData:
            print("CorrelationData")
            results[path][k] = {}
            results[path][k]['path'] = v.path
            results[path][k]['type'] = 'heatmap'
        else:
            pass

    return results
