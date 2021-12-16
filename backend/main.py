from fastapi import FastAPI, Response
import uvicorn
from starlette.middleware.cors import CORSMiddleware
import sys
sys.path.append('../optinist')
from routers import files, run, params, outputs, algolist

app = FastAPI(docs_url="/docs", openapi_url="/openapi")
app.include_router(algolist.router)
app.include_router(files.router)
app.include_router(outputs.router)
app.include_router(params.router)
app.include_router(run.router)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"]
)

@app.get("/")
async def root():
    return {"message": "Hello World"}

<<<<<<< HEAD
@app.get("/api/params/{name}")
async def params(name: str):
    config = {}
    filepath = f'../optinist/config/{name}.yaml'
    if os.path.exists(filepath):
        with open(filepath) as f:
            config = yaml.safe_load(f)
    print(config)
    return config

def get_nest_dict(value):
    algo_dict = {}
    for _k, _v in value.items():
        algo_dict[_k] = {}
        if type(_v) is dict:
            algo_dict[_k]['children'] = get_nest_dict(_v)
        else:
            algo_dict[_k] = {}

            # get args
            sig = inspect.signature(_v)
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

    return algo_dict


@app.get("/api/algolist")
async def run() -> List:
    {
        'caiman_mc': {
            'args': ['images', 'timeseries']
        },
        'caiman_cnmf': {
            'args': ['images', 'timeseries']
        }
    }

    algo_dict = get_nest_dict(wrapper_dict)

    return algo_dict

@app.get("/api/cookie-test")
=======
@app.get("/cookie-test")
>>>>>>> f8997191e038285677ae729c9bd7b33d2b12bda4
def create_cookie(response: Response):
    response.set_cookie(key="fakesession", value="fake-cookie-session-value")
    return {"message": "cookie is set."}

if __name__ == '__main__':
	uvicorn.run('main:app', host='0.0.0.0', port=8000, reload=True)
