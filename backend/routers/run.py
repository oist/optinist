import traceback
import os
import gc
import json
import tracemalloc
import time
# from pytools.persistent_dict import PersistentDict

from wrappers.data_wrapper import *
from wrappers.optinist_exception import AlgorithmException

from .const import BASE_DIR, OPTINIST_DIR
from .utils.save import save_nwb
from .utils.run import run_algorithm
from .utils.memory import display_top
from .utils.results import get_results
from .utils.utils import algo_network
from .utils.set_data import set_data
from .utils.snakemake import create_snakemake_files

from fastapi import APIRouter, WebSocket
router = APIRouter()


@router.get("/run/ready/{request_id}")
async  def run_ready(request_id: str):
    return {'message': 'ready...', 'status': 'ready', "requestId": request_id}


@router.websocket("/run")
async def websocket_endpoint(websocket: WebSocket):
    # import run_pipeline

    try:
        await websocket.accept()
        # Wait for any message from the client
        message = await websocket.receive_text()
        graph, startNodeList, nodeDict, edgeList, endNodeList = algo_network(message)

        tracemalloc.start()

        ### snakemake
        all_outputs = create_snakemake_files(BASE_DIR, OPTINIST_DIR, nodeDict, edgeList, endNodeList)

        error_flag = False
        time.sleep(1)
        # Send message to the client
        for key, value in all_outputs.items():
            info = value["info"]
            path = value["path"]
            label = value['label']

            if isinstance(info, str) or isinstance(info, list):
                print("error")
                if isinstance(info, str):
                    error_message = info
                else:
                    error_message = "¥n".join(info)
                error_label = label
                error_flag = True
            else:
                print(f"finish {label}")
                results = get_results(info, path)
                time.sleep(5)
                await websocket.send_json({
                    'message': label + ' success',
                    'status': 'success',
                    'outputPaths': results
                })
                del results

        gc.collect()

        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics("lineno")
        display_top(snapshot, top_stats)

        if error_flag:
            time.sleep(5)
            await websocket.send_json({
                'message': error_message,
                'name': error_label,
                'status': 'error'
            })
        else:
            time.sleep(5)
            await websocket.send_json({
                'message': 'completed',
                'status': 'completed'
            })

    except AlgorithmException as e:
        time.sleep(1)
        await websocket.send_json({
            'message': e.get_message(),
            # 'name': item['data']['label'],
            'status': 'error'
        })
        await websocket.close()
    except Exception as e:
        message  = list(traceback.TracebackException.from_exception(e).format())[-1]
        print(traceback.format_exc())
        time.sleep(1)
        await websocket.send_json({
            'message': message,
            # 'name': item['data']['label'],
            'status': 'error'
        })
        await websocket.close()
    finally:
        gc.collect()
        print('Bye..')
        await websocket.close()


def add_other_input(graph, item, nodeDict, info):
    other_input = None
    if len(graph[item['id']]) > 0:
        prev_id = item["id"]
        node_id = list(graph[item['id']].keys())[0]
        for k, v in graph.items():
            if k != prev_id and k != node_id and node_id in v:
                item = nodeDict[k]
                item['key'] = list(graph[k].values())[0]
                other_input = set_data(None, item, None)

    if other_input is not None:
        for k, v in other_input.items():
            if k in info.keys():
                info[f"{k}{np.random.randint(10000)}"] = v
            else:
                info[k] = v
    
    return info
