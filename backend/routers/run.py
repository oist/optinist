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
    flowList = list(map(lambda x: FlowItem(**x), json.loads(flowList)))
    try:
        for item in flowList:

            await websocket.send_json({'message': item.label+' started', 'status': 'ready'})

            # run algorithm
            info = None
            if item.type == 'image':
                info = {'path': ImageData(item.path, '')}
            elif item.type == 'algo':
                # parameterをint, floatに変換
                from .utils import string_to_float
                item.param = string_to_float(item.param)

                wrapper = get_dict_leaf_value(wrapper_dict, item.path.split('/'))
                info = wrapper["function"](
                    *prev_info.values(), params=item.param)

            prev_info = info

            assert info is not None

            results = OrderedDict()
            results[item.path] = {}
            for k, v in info.items():
                if type(v) is ImageData:
                    print("ImageData")
                    results[item.path][k] = {}
                    results[item.path][k]['path'] = v.json_path
                    results[item.path][k]['type'] = 'images'
                    results[item.path][k]['max_index'] = len(v.data)
                elif type(v) is TimeSeriesData:
                    print("TimeSeriesData")
                    results[item.path][k] = {}
                    results[item.path][k]['path'] = v.path
                    results[item.path][k]['type'] = 'timeseries'
                elif type(v) is CorrelationData:
                    print("CorrelationData")
                    results[item.path][k] = {}
                    results[item.path][k]['path'] = v.path
                    results[item.path][k]['type'] = 'heatmap'
                else:
                    pass

            # Send message to the client
            await websocket.send_json({'message': item.label+' success', 'status': 'success', 'outputPaths': results})

        await websocket.send_json({'message': 'completed', 'status': 'completed'})
    except AlgorithmException as e:
        await websocket.send_json({'message': e.get_message(), 'name': item.label, 'status': 'error'})
    except Exception as e:
        message  = list(traceback.TracebackException.from_exception(e).format())[-1]
        print(traceback.format_exc())
        await websocket.send_json({'message': message, 'name': item.label, 'status': 'error'})
    finally:
        print('Bye..')
        await websocket.close()
