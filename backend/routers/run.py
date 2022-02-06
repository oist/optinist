
from fastapi import APIRouter, BackgroundTasks
from pydantic import BaseModel
from typing import List

from workflow.params import get_snakemake_params
from cui_api.snakemake import run_snakemake


router = APIRouter()


class RunItem(BaseModel):
    nodeList: list = []
    edgeList: list = []
    snakemakeParam: dict = {}
    nwbParam: dict = {}


@router.post("/run")
async def params(runItem: RunItem, background_tasks: BackgroundTasks):
    import uuid
    unique_id = str(uuid.uuid4())

    from workflow.set_pipeline import set_pipeline
    set_pipeline(runItem)
    
    snakemake_params = get_snakemake_params(runItem.snakemakeParam)
    background_tasks.add_task(run_snakemake, snakemake_params)

    # timestamp = "20220205"
    # background_tasks.add_task(dummy_run_pipeline, unique_id)
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
            output[info["nodeId"]] = {
                "status": info["status"],
                "message": info["message"],
                "name": info["name"],
                "outputPaths": info["outputPaths"],
            }
        else:
            output[info["nodeId"]] = {
                "status": info["status"],
                "message": info["message"],
                "name": info["name"],
            }

    return output


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
