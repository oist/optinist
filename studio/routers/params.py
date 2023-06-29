from typing import Any, Dict

from fastapi import APIRouter

from studio.core.utils.config_handler import ConfigReader
from studio.core.utils.filepath_finder import find_param_filepath
from studio.schemas.params import SnakemakeParams

router = APIRouter()


@router.get("/params/{name}", response_model=Dict[str, Any], tags=["params"])
async def get_params(name: str):
    filepath = find_param_filepath(name)
    config = ConfigReader.read(filepath)
    return config


@router.get("/snakemake", response_model=SnakemakeParams, tags=["params"])
async def get_snakemake_params():
    filepath = find_param_filepath("snakemake")
    return ConfigReader.read(filepath)
