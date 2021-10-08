import os

import imageio
from PIL import Image
from fastapi import Depends, FastAPI, File, Response, UploadFile
from fastapi.staticfiles import StaticFiles
import uvicorn
from starlette.middleware.cors import CORSMiddleware
from typing import List, Optional
from pydantic import BaseModel
import sys
import yaml
sys.path.append('../optinist')
from wrappers import wrapper_dict

app = FastAPI()

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

@app.get("/")
async def root():
    return {"message": "Hello World"}

@app.get("/params/{name}")
async def params(name: str):
    with open(f'../optinist/config/{name}.yaml') as f:
        config = yaml.safe_load(f)
    print(config)
    return config

@app.get("/algolist")
async def run() -> List:
    print(wrapper_dict.keys())
    return list(wrapper_dict.keys())

@app.get("/cookie-test")
def create_cookie(response: Response):
    response.set_cookie(key="fakesession", value="fake-cookie-session-value")
    return {"message": "cookie is set."}

app.mount("/files", StaticFiles(directory="files"), name="files")

@app.post("/upload/{path}")
async def create_file(response: Response, path: str, file: UploadFile = File(...)):
    max_index = 30
    folder_name = path
    contents = await file.read()
    file.filename = "tmp.tiff"
    with open(os.path.join("_tmp", file.filename), "wb") as f:
        f.write(contents)

    os.makedirs(os.path.join("files", folder_name), exist_ok=True)
    tiffs = imageio.volread(os.path.join("_tmp", file.filename))
    counter = 0
    for i, tiff_data in enumerate(tiffs):
        if counter == max_index:
            break
        img = Image.fromarray(tiff_data)
        img = img.convert("L")
        img.save(os.path.join("files", folder_name, f"{i}.png"))
        counter += 1

    os.remove(os.path.join("_tmp", file.filename))

    response.set_cookie(key="directory", value=folder_name)

    return {"folderName": folder_name, "maxIndex": max_index}


@app.post("/run")
async def run(flowList: List[FlowItem]):
    print('run_code')
    print(wrapper_dict)
    import run
    return run.run_code(wrapper_dict, flowList)


if __name__ == '__main__':
	uvicorn.run('main:app', host='0.0.0.0', port=8000, reload=True)
