import os

import imageio
from PIL import Image
from fastapi import Depends, FastAPI, File, Response, UploadFile, Form
from fastapi.staticfiles import StaticFiles
import uvicorn
from starlette.middleware.cors import CORSMiddleware
from typing import Dict, List, Optional
from pydantic import BaseModel
import sys
import yaml
import json
sys.path.append('../optinist')
from wrappers import wrapper_dict

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
    param: Optional[Dict] = None

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

@app.get("/api/algolist")
async def run() -> List:
    print(wrapper_dict.keys())
    return list(wrapper_dict.keys())

@app.get("/api/cookie-test")
def create_cookie(response: Response):
    response.set_cookie(key="fakesession", value="fake-cookie-session-value")
    return {"message": "cookie is set."}

os.makedirs('files', exist_ok=True)
app.mount("/api/files", StaticFiles(directory="files"), name="files")

@app.post("/api/upload/{fileName}/{inputFileNumer}")
async def create_file(response: Response, fileName: str, element_id: str = Form(...), file: UploadFile = File(...), inputFileNumer: int=1):
    # max_index = 30
    root_folder = os.path.join("files", fileName+"("+element_id+")")
    png_folder = os.path.join(root_folder, "pngs")
    tiff_folder = os.path.join(root_folder, "tiff")
    os.makedirs(root_folder, exist_ok=True)
    os.makedirs(png_folder, exist_ok=True)
    os.makedirs(tiff_folder, exist_ok=True)
    contents = await file.read()
    file.filename = fileName
    tiff_path = os.path.join(tiff_folder, file.filename)
    with open(tiff_path, "wb") as f:
        f.write(contents)

    tiffs = imageio.volread(tiff_path)[:inputFileNumer]

    for i, tiff_data in enumerate(tiffs):
        img = Image.fromarray(tiff_data)
        img = img.convert("L")
        img.save(os.path.join(png_folder, f"{i}.png"))

    response.set_cookie(key="directory", value=png_folder)

    return {"pngFolder": png_folder, "tiffPath": tiff_path, "maxIndex": len(tiffs)}


@app.post("/api/run")
async def run(flowList: List[FlowItem]):
    print('run_code')
    print(wrapper_dict)
    import run
    return run.run_code(wrapper_dict, flowList)

@app.get("/api/outputs/{file_path:path}")
async def read_file(file_path: str):
    with open(os.path.join(".", file_path), 'r') as f:
        json_dict = json.load(f)
    return { "data": json_dict }

if __name__ == '__main__':
	uvicorn.run('main:app', host='0.0.0.0', port=8000, reload=True)
