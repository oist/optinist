
from fastapi import APIRouter
import os
import yaml
from .utils.utils import dict2nest
from .utils.params import get_params

router = APIRouter()


@router.get("/nwb")
async def params():
    filepath = os.path.join('..', 'optinist', 'config', 'nwb.yaml')
    config = get_params(filepath)
    
    config = dict2nest(config)

    return config
