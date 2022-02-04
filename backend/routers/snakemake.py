
from fastapi import APIRouter
import os
from .utils.params import get_params

router = APIRouter()


@router.get("/snakemake")
async def params():
    filepath = os.path.join('..', 'optinist', 'config', 'snakemake.yaml')
    config = get_params(filepath)
    return config
