from fastapi import Depends, FastAPI
import uvicorn
from starlette.middleware.cors import CORSMiddleware
from typing import List, Optional
from pydantic import BaseModel
import sys
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
    if name == 'test':
        import params
        return params.get_params()

@app.get("/algolist")
async def run() -> List:
    print(wrapper_dict.keys())
    return list(wrapper_dict.keys())

@app.post("/run")
async def run(flowList: List[FlowItem]):
    print('run_code')
    print(wrapper_dict)
    import run
    return run.run_code(wrapper_dict, flowList)


if __name__ == '__main__':
	uvicorn.run('main:app', host='0.0.0.0', port=8000, reload=True)
