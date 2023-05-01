from typing import Any, Dict

from fastapi import APIRouter

from optinist.api.config.config_reader import ConfigReader
from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import join_filepath
from optinist.routers.model import NWBParams, SnakemakeParams

router = APIRouter()


@router.get("/params/{name}", response_model=Dict[str, Any], tags=["params"])
async def get_params(name: str):
    filepath = join_filepath([DIRPATH.CONFIG_DIR, f"{name}.yaml"])
    config = ConfigReader.read(filepath)
    return config


@router.get("/snakemake", response_model=SnakemakeParams, tags=["params"])
async def get_snakemake_params():
    filepath = join_filepath([DIRPATH.CONFIG_DIR, "snakemake.yaml"])
    return ConfigReader.read(filepath)


@router.get("/nwb", response_model=NWBParams, tags=["params"])
async def get_nwb_params():
    filepath = join_filepath([DIRPATH.CONFIG_DIR, "nwb.yaml"])
    return ConfigReader.read(filepath)
