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
import pickle


router = APIRouter()


@router.get("/run/ready/{request_id}")
async  def run_ready(request_id: str):
    return {'message': 'ready...', 'status': 'ready', "requestId": request_id}



class RunItem(BaseModel):
    nodeList: list = []
    edgeList: list = []
    snakemakeParam: dict = {}
    nwbParam: dict = {}


def dummy_run_pipeline(unique_id):
    import time
    i = 0
    os.makedirs(f"/tmp/optinist/{unique_id}")
    while True:
        print(f"i = {str(i)}")
        if i > 60:
            break

        if i == 5:
            with open(f"/tmp/optinist/{unique_id}/A.pkl", "wb") as f:
                info = {
                    'status': 'success',
                    'message': 'A success',
                    'name': 'A',
                    'outputPaths': {
                        "A_image": {
                            "path": f"/tmp/optinist/{unique_id}/A_images.json",
                            "type": "images",
                        },
                        "A_timeseries": {
                            "path": f"/tmp/optinist/{unique_id}/A_timeseries.json",
                            "type": "timeseries",
                        },
                    }
                }
                pickle.dump(info, f)
            

        if i == 10:
            with open(f"/tmp/optinist/{unique_id}/B.pkl", "wb") as f:
                info = {
                    'status': 'error',
                    'message': 'error reason',
                    'name': 'B',
                }
                pickle.dump(info, f)

        if i == 15:
            with open(f"/tmp/optinist/{unique_id}/C.pkl", "wb") as f:
                info = {
                    'status': 'success',
                    'message': 'C success',
                    'name': 'C',
                    'outputPaths': {
                        "C_image": {
                            "path": f"/tmp/optinist/{unique_id}/C_heatmp.json",
                            "type": "heatmap",
                        },
                        "C_timeseries": {
                            "path": f"/tmp/optinist/{unique_id}/C_timeseries.json",
                            "type": "timeseries",
                        },
                    }
                }
                pickle.dump(info, f)

        i += 1
        time.sleep(1)


@router.post("/run")
async def params(runItem: RunItem, background_tasks: BackgroundTasks):
    import uuid
    unique_id = str(uuid.uuid4())
    timestamp = "20220205"
    background_tasks.add_task(dummy_run_pipeline, unique_id)
    print("start hevy task")

    return unique_id


@router.post("/run/result/{uid}")
async def params(uid: str):
    from glob import glob
    print(f"/tmp/optinist/{uid}/*.pkl")
    runPaths = glob(f"/tmp/optinist/{uid}/*.pkl")
    print(runPaths)

    output = {}
    for request_path in runPaths:
        with open(request_path, "rb") as f:
            info = pickle.load(f)
        print(info)

        if info["status"] == "success":
            output[request_path] = {
                "status": info["status"],
                "message": info["message"],
                "name": info["name"],
                "outputPaths": info["outputPaths"],
            }
        else:
            output[request_path] = {
                "status": info["status"],
                "message": info["message"],
                "name": info["name"],
            }

    return output


# @router.post("/run/path")
# async def params(runPaths: List[str] = []):
#     run_success = {
#         'savePaths': ["/tmp/optinist/0/A.out"],
#         'status': 'success',
#         'outputPaths': ["A.json"]
#     }
#     return run_success


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
