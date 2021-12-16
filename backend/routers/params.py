
from fastapi import APIRouter
import os
import yaml

router = APIRouter()

@router.get("/params/{name}")
async def params(name: str):
    config = {}
    filepath = os.path.join('..', 'optinist', 'config', f'{name}.yaml')
    print('filepath:', filepath)
    if os.path.exists(filepath):
        with open(filepath) as f:
            config = yaml.safe_load(f)
    print(config)
    return config
