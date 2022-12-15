from fastapi import APIRouter

from studio.api.config.config_reader import ConfigReader
from studio.api.dir_path import DIRPATH
from studio.api.utils.filepath_creater import join_filepath

router = APIRouter()


@router.get("/params/{name}")
async def params(name: str):
    filepath = join_filepath([DIRPATH.CONFIG_DIR, f'{name}.yaml'])
    config = ConfigReader.read(filepath)
    return config


@router.get("/snakemake")
async def snakemake_params():
    filepath = join_filepath([DIRPATH.CONFIG_DIR, f'snakemake.yaml'])
    return ConfigReader.read(filepath)


@router.get("/nwb")
async def nwb_params():
    filepath = join_filepath([DIRPATH.CONFIG_DIR, f'nwb.yaml'])
    return ConfigReader.read(filepath)
