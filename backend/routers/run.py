import traceback
import os
import gc
import json
import tracemalloc
import time

from wrappers.data_wrapper import *
from wrappers.optinist_exception import AlgorithmException

from .const import BASE_DIR, OPTINIST_DIR
from .utils.memory import display_top
from .utils.results import get_results
from .utils.utils import algo_network
from .utils.snakemake import create_snakemake_files

from fastapi import APIRouter, BackgroundTasks #, WebSocket
from pydantic import BaseModel
import sys
from typing import List


router = APIRouter()


@router.get("/run/ready/{request_id}")
async  def run_ready(request_id: str):
    return {'message': 'ready...', 'status': 'ready', "requestId": request_id}



class RunItem(BaseModel):
    uid: str
    timestamp: str
    nodeList: list = []
    edgeList: list = []
    snakemakeParam: dict = {}
    nwbParam: dict = {}


savePaths = {}

def dummy_run_pipeline():
    import time
    i = 0
    while True:
        print(f"i = {str(i)}")
        if i > 60:
            break

        if i == 5:
            savePaths["/tmp/optinist/0/A.out"] = {
                "path": "A.json",
                "status": "success"
            }

        if i == 10:
            savePaths["/tmp/optinist/0/B.out"] = {
                "path": None,
                "status": "error"
            }

        if i == 15:
            savePaths["/tmp/optinist/0/C.out"] = {
                "path": "C.json",
                "status": "success"
            }

        i+=1
        time.sleep(1)


@router.post("/run")
async def params(runItem: RunItem, background_tasks: BackgroundTasks):
    runPaths = ["/tmp/optinist/0/A.out", "/tmp/optinist/0/B.out", "/tmp/optinist/0/C.out"]
    background_tasks.add_task(dummy_run_pipeline)
    print("finish hevy task")

    return runPaths


@router.post("/run/path")
async def params(runPaths: List[str] = []):
    """
        output = {
            /tmp/optinist/0/A.out: {
                "path": A.json
                "status": "success"
            },
            /tmp/optinist/0/B.out: {
                "path": None,
                "status": "error"
            },
            /tmp/optinist/0/C.out: {
                "path": "C.json",
                "status": "success"
            }
        }
    """
    output = {}
    print(savePaths)
    for request_path in runPaths:
        if request_path in savePaths:
            output[request_path] = {
                "path": savePaths[request_path]["path"],
                "status": savePaths[request_path]["status"]
            }

    return output


@router.post("/run/path")
async def params(runPaths: List[str] = []):
    run_success = {
        'savePaths': ["/tmp/optinist/0/A.out"],
        'status': 'success',
        'outputPaths': ["A.json"]
    }
    return run_success


# @router.websocket("/run")
# async def websocket_endpoint(websocket: WebSocket):
#     # import run_pipeline

#     try:
#         await websocket.accept()
#         # Wait for any message from the client
#         message = await websocket.receive_text()
#         message = json.loads(message)

#         tracemalloc.start()

#         ### snakemake
#         all_outputs = create_snakemake_files(BASE_DIR, OPTINIST_DIR, message)

#         error_flag = False
#         time.sleep(1)
#         # Send message to the client
#         for key, value in all_outputs.items():
#             info = value["info"]
#             path = value["path"]
#             label = value['label']

#             if isinstance(info, str) or isinstance(info, list):
#                 print("error")
#                 if isinstance(info, str):
#                     error_message = info
#                 else:
#                     error_message = "¥n".join(info)
#                 error_label = label
#                 error_flag = True
#             else:
#                 print(f"finish {label}")
#                 results = get_results(info, path)
#                 time.sleep(1)
#                 await websocket.send_json({
#                     'message': label + ' success',
#                     'status': 'success',
#                     'outputPaths': results
#                 })
#                 del results

#         gc.collect()

#         snapshot = tracemalloc.take_snapshot()
#         top_stats = snapshot.statistics("lineno")
#         # display_top(snapshot, top_stats)

#         if error_flag:
#             time.sleep(1)
#             await websocket.send_json({
#                 'message': error_message,
#                 'name': error_label,
#                 'status': 'error'
#             })
#         else:
#             time.sleep(1)
#             await websocket.send_json({
#                 'message': 'completed',
#                 'status': 'completed'
#             })

#     except AlgorithmException as e:
#         time.sleep(1)
#         await websocket.send_json({
#             'message': e.get_message(),
#             # 'name': item['data']['label'],
#             'status': 'error'
#         })
#         await websocket.close()
#     except Exception as e:
#         message  = list(traceback.TracebackException.from_exception(e).format())[-1]
#         # print(traceback.format_exc())
#         time.sleep(1)
#         await websocket.send_json({
#             'message': message,
#             # 'name': item['data']['label'],
#             'status': 'error'
#         })
#         await websocket.close()
#     finally:
#         gc.collect()
#         print('Bye..')
#         await websocket.close()
