import os
import inspect

import imageio
from PIL import Image
from fastapi import Depends, FastAPI, File, Response, UploadFile, Form, WebSocket
from fastapi.staticfiles import StaticFiles
import uvicorn
from starlette.middleware.cors import CORSMiddleware
from typing import Dict, List, Optional
from pydantic import BaseModel
import sys
import yaml
import json
import pandas as pd
sys.path.append('../optinist')
from wrappers import wrapper_dict
from collections import OrderedDict
from wrappers.data_wrapper import *
from wrappers.optinist_exception import AlgorithmException
import traceback

app = FastAPI(docs_url="/api/docs", openapi_url="/api")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"]
)

class FlowItem(BaseModel):
    label: str
    path: Optional[str] = None
    type: str
    param: Optional[Dict] = {}

@app.get("/")
async def root():
    return {"message": "Hello World"}

@app.get("/api/params/{name}")
async def params(name: str):
    config = {}
    filepath = f'../optinist/config/{name}.yaml'
    if os.path.exists(filepath):
        with open(filepath) as f:
            config = yaml.safe_load(f)
    print(config)
    return config

def get_nest_dict(value, parent_k):
    algo_dict = {}
    for _k, _v in value.items():
        algo_dict[_k] = {}
        if type(_v) is dict and 'function' not in _v.keys():
            algo_dict[_k]['children'] = get_nest_dict(
                _v, parent_k+'/'+_k if parent_k != '' else _k)
        else:
            algo_dict[_k] = {}

            # get args
            sig = inspect.signature(_v['function'])
            algo_dict[_k]['args'] = [
                {
                    'name': x.name, 
                    'type': x.annotation.__name__
                }
                for x in sig.parameters.values()
                if x.name != 'params'
            ]

            # get returns
            if sig.return_annotation is not inspect._empty:
                algo_dict[_k]['returns'] = [
                    {
                        'name': k,
                        'type': v.__name__
                    }
                    for k, v in sig.return_annotation.items()
                ]
            
            # parameter path
            if 'parameter' in _v.keys():
                algo_dict[_k]['parameter'] = _v['parameter']
            else:
                algo_dict[_k]['parameter'] = ''
            
            # path
            algo_dict[_k]['path'] = parent_k + '/' + _k

    return algo_dict

def get_dict_leaf_value(root_dict: dict, path_list: List[str]):
    path = path_list.pop(0)
    if(len(path_list) > 0):
        return get_dict_leaf_value(root_dict[path],path_list)
    else:
        return root_dict[path]

@app.get("/api/algolist")
async def run() -> List:
    # print(wrapper_dict.keys())
    {
        'caiman': {
            'children': {
                'caiman_mc' : {
                    'args': ['images', 'timeseries'],
                    'path': 'caiman/caiman_mc'
                },
                'caiman_cnmf': {
                    'args': ['images', 'timeseries'],
                    'path': 'caiman/caiman_mc'
                }
            }
        }
    }
    {
        'caiman.caiman_mc': {

        },
        'caiman.caiman_cnmf': {

        },
    }

    algo_dict = get_nest_dict(wrapper_dict, '')

    return algo_dict

@app.get("/api/cookie-test")
def create_cookie(response: Response):
    response.set_cookie(key="fakesession", value="fake-cookie-session-value")
    return {"message": "cookie is set."}

os.makedirs('files', exist_ok=True)
app.mount("/api/files", StaticFiles(directory="files"), name="files")

@app.post("/api/upload/{fileName}/{inputFileNumer}")
async def create_file(response: Response, fileName: str, element_id: str = Form(...), file: UploadFile = File(...), inputFileNumer: int=1):
    root_dir = os.path.join("files", fileName+"("+element_id+")")
    os.makedirs(root_dir, exist_ok=True)

    contents = await file.read()
    tiff_file_path = os.path.join(root_dir, fileName)
    with open(tiff_file_path, "wb") as f:
        f.write(contents)

    tiffs = imageio.volread(tiff_file_path)[:inputFileNumer]

    images = []
    for i, _img in enumerate(tiffs):
        images.append(_img.tolist())

    json_data_path = os.path.join(root_dir, 'image.json')
    pd.DataFrame(images).to_json(json_data_path, indent=4, orient="values")

    return {"json_data_path": json_data_path, "tiff_file_path": tiff_file_path}

@app.get("/api/run/ready/{request_id}")
async  def run_ready(request_id: str):
    return {'message': 'ready...', 'status': 'ready', "requestId": request_id}

@app.websocket("/run")
async def websocket_endpoint(websocket: WebSocket):
    import json
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
            if item.type == 'data':
                info = {'path': ImageData(item.path, '')}
            elif item.type == 'algo':
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
                    results[item.path][k]['path'] = v.path
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

@app.get("/api/outputs/{file_path:path}")
async def read_file(file_path: str):
    with open(os.path.join(".", file_path), 'r') as f:
        json_dict = json.load(f)
    return { "data": json_dict }

if __name__ == '__main__':
	uvicorn.run('main:app', host='0.0.0.0', port=8000, reload=True)
