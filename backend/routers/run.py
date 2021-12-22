from fastapi import APIRouter, WebSocket
from pydantic import BaseModel
from typing import Dict, List, Optional
import traceback
import sys

from wrappers import wrapper_dict
from collections import OrderedDict
from wrappers.data_wrapper import *
from wrappers.optinist_exception import AlgorithmException
from .utils import get_algo_network, run_algorithm, get_results

sys.path.append('../../optinist')

router = APIRouter()



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
