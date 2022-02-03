
from fastapi import APIRouter
import os
import yaml
from .utils.params import get_params

router = APIRouter()


@router.get("/nwb")
async def params():
    filepath = os.path.join('..', 'optinist', 'config', 'nwb.yaml')
    config = get_params(filepath)

    return config
