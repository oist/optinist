from typing import Any, Dict

from fastapi import APIRouter

from studio.app.common.core.workflow.workflow_params import read_default_params
from studio.app.common.schemas.params import SnakemakeParams

router = APIRouter(tags=["params"])


@router.get("/params/{name}", response_model=Dict[str, Any])
async def get_params(name: str):
    return read_default_params(name)


@router.get("/snakemake", response_model=SnakemakeParams)
async def get_snakemake_params():
    return read_default_params("snakemake")
