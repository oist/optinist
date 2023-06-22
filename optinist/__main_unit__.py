import argparse
import os

import uvicorn
from fastapi import FastAPI, Request
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from starlette.middleware.cors import CORSMiddleware

from optinist.api.config.config_reader import ConfigReader
from optinist.api.dir_path import DIRPATH as OPTINIST_DIRPATH
from optinist.api.utils.filepath_creater import join_filepath
from optinist.routers import algolist, experiment, files, hdf5, outputs, params, run

app = FastAPI(docs_url="/docs", openapi_url="/openapi")
app.include_router(algolist.router)
app.include_router(files.router)
app.include_router(outputs.router)
app.include_router(params.router)
app.include_router(run.router)
app.include_router(hdf5.router)
app.include_router(experiment.router)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

if "OPTINIST_DEV_ROOT_DIR" in os.environ:
    FRONTEND_DIRPATH = os.environ["OPTINIST_DEV_ROOT_DIR"] + "/frontend"
else:
    FRONTEND_DIRPATH = OPTINIST_DIRPATH.ROOT_DIR + "/frontend"

app.mount(
    "/static",
    StaticFiles(directory=f"{FRONTEND_DIRPATH}/build/static"),
    name="static",
)

templates = Jinja2Templates(directory=f"{FRONTEND_DIRPATH}/build")


@app.get("/")
async def root(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})


def main(develop_mode: bool = False):
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", type=str, default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--reload", action="store_true")
    args = parser.parse_args()

    # set fastapi@uvicorn logging config.
    logging_config_path = join_filepath(
        [OPTINIST_DIRPATH.ROOT_DIR, "app_config", "logging.yaml"]
    )
    logging_config = ConfigReader.read(logging_config_path)
    fastapi_logging_config = uvicorn.config.LOGGING_CONFIG
    fastapi_logging_config["formatters"]["default"]["fmt"] = logging_config[
        "fastapi_logging_config"
    ]["default_fmt"]
    fastapi_logging_config["formatters"]["access"]["fmt"] = logging_config[
        "fastapi_logging_config"
    ]["access_fmt"]

    if develop_mode:
        reload_options = {"reload_dirs": ["optinist"]} if args.reload else {}
        uvicorn.run(
            "optinist.__main_unit__:app",
            host=args.host,
            port=args.port,
            log_config=fastapi_logging_config,
            reload=args.reload,
            **reload_options,
        )
    else:
        uvicorn.run(
            "optinist.__main_unit__:app",
            host=args.host,
            port=args.port,
            log_config=fastapi_logging_config,
            reload=False,
        )
