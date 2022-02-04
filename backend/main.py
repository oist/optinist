from fastapi import FastAPI, Response
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
    snakemake
)

app = FastAPI(docs_url="/docs", openapi_url="/openapi")
app.include_router(algolist.router)
app.include_router(files.router)
app.include_router(outputs.router)
app.include_router(params.router)
app.include_router(run.router)
app.include_router(nwb.router)
app.include_router(snakemake.router)

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

@app.get("/cookie-test")
def create_cookie(response: Response):
    response.set_cookie(key="fakesession", value="fake-cookie-session-value")
    return {"message": "cookie is set."}

if __name__ == '__main__':
	uvicorn.run('main:app', host='0.0.0.0', port=8000, reload=True)
