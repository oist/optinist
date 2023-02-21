from typing import Dict, Any
from fastapi import APIRouter

from optinist.api.config.config_reader import ConfigReader
from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import join_filepath
from optinist.routers.model import NWB as  NWBModel, Snakemake as SnakemakeModel

router = APIRouter()


@router.get("/params/{name}", response_model=Dict[str, Any])
async def params(name: str):
    filepath = join_filepath([DIRPATH.CONFIG_DIR, f'{name}.yaml'])
    config = ConfigReader.read(filepath)
    return config


@router.get("/snakemake", response_model=SnakemakeModel)
async def snakemake_params():
    filepath = join_filepath([DIRPATH.CONFIG_DIR, f'snakemake.yaml'])
    return ConfigReader.read(filepath)


@router.get("/nwb", response_model=NWBModel)
async def nwb_params():
    filepath = join_filepath([DIRPATH.CONFIG_DIR, f'nwb.yaml'])
    return ConfigReader.read(filepath)
