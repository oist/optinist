import traceback
import os
import gc
import json
import tracemalloc

from wrappers.data_wrapper import *
from wrappers.optinist_exception import AlgorithmException

from .const import BASE_DIR
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

    await websocket.accept()
    # Wait for any message from the client

    message = await websocket.receive_text()
    graph, startNodeList, nodeDict = algo_network(message)
    create_snakemake_files(nodeDict)

    message = json.loads(message)

    nwb_params = message['nwbParam']

    try:
        tracemalloc.start()
        item = nodeDict[startNodeList[0]]
        item['key'] = list(graph[startNodeList[0]].values())[0]
        prev_info = None
        info = None
        results = None
        node_id = None

        while True:

            await websocket.send_json({
                'message': item['data']['label'] + ' started',
                'status': 'ready'
            })

            # 実行
            if item['type'] == 'AlgorithmNode':
                info = run_algorithm(info, item)
            else:
                info = set_data(info, item, nwb_params)

            # prev_infoを登録
            if 'nwbfile' not in info.keys():
                info['nwbfile'] = prev_info['nwbfile']
            else:
                assert info is not None and 'nwbfile' in info.keys()

            results = get_results(info, item)

            # Send message to the client
            await websocket.send_json({
                'message': item['data']['label'] + ' success',
                'status': 'success',
                'outputPaths': results
            })

            """ add other input """
            info = add_other_input(graph, item, nodeDict, info)
            """ add other input """

            if len(graph[item['id']]) > 0:
                node_id = list(graph[item['id']].keys())[0]
                item = nodeDict[node_id]
            else:
                break

            prev_info = info

        # NWBを保存
        if info['nwbfile'] is not None:
            save_nwb(info['nwbfile'], os.path.join(BASE_DIR, 'nwb'))

        del prev_info, info, results
        gc.collect()

        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics("lineno")
        display_top(snapshot, top_stats)

        await websocket.send_json({
            'message': 'completed',
            'status': 'completed'
        })

    except AlgorithmException as e:
        await websocket.send_json({
            'message': e.get_message(),
            'name': item['data']['label'],
            'status': 'error'
        })
    except Exception as e:
        message  = list(traceback.TracebackException.from_exception(e).format())[-1]
        print(traceback.format_exc())
        await websocket.send_json({
            'message': message,
            'name': item['data']['label'],
            'status': 'error'
        })
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
