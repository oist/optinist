
from fastapi import APIRouter
import os
import yaml
from .utils.utils import dict2nest
from .utils.params import get_params

router = APIRouter()


@router.get("/snakemake")
async def params():
    filepath = os.path.join('..', 'optinist', 'config', 'snakemake.yaml')
    config = get_params(filepath)
    
    config = dict2nest(config)

    return config
