from typing import Any, Dict

from fastapi import APIRouter

from studio.app.common.core.utils.config_handler import ConfigReader
from studio.app.common.core.utils.filepath_finder import find_param_filepath
from studio.app.common.schemas.params import SnakemakeParams

router = APIRouter(tags=["params"])


@router.get("/params/{name}", response_model=Dict[str, Any])
async def get_params(name: str):
    filepath = find_param_filepath(name)
    config = ConfigReader.read(filepath)
    return config


@router.get("/snakemake", response_model=SnakemakeParams)
async def get_snakemake_params():
    filepath = find_param_filepath("snakemake")
    return ConfigReader.read(filepath)
