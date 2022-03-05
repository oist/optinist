from fastapi import FastAPI, Response, Request
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

import uvicorn
from starlette.middleware.cors import CORSMiddleware
import sys
sys.path.append('../optinist')
from routers import (
    files,
    run,
    params,
    outputs,
    algolist,
    nwb,
    snakemake,
    hdf5,
    experiment,
)


app = FastAPI(docs_url="/docs", openapi_url="/openapi")
app.include_router(algolist.router)
app.include_router(files.router)
app.include_router(outputs.router)
app.include_router(params.router)
app.include_router(run.router)
app.include_router(nwb.router)
app.include_router(snakemake.router)
app.include_router(hdf5.router)
app.include_router(experiment.router)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"]
)

app.mount(
    "/static",
    StaticFiles(directory="../frontend/build/static"),
    name="static"
)

templates = Jinja2Templates(directory="../frontend/build")

@app.get("/")
async def root(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})

if __name__ == '__main__':
    uvicorn.run('main:app', host='0.0.0.0', port=8000, reload=True)
