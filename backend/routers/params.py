import os
from fastapi import APIRouter
from .utils.params import get_params

router = APIRouter()


@router.get("/params/{name}")
async def params(name: str):
    config = {}
    filepath = os.path.join('..', 'optinist', 'config', f'{name}.yaml')
    config = get_params(filepath)
    return config
