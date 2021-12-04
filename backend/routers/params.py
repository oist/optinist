
from fastapi import APIRouter
import os
import yaml

router = APIRouter()

@router.get("/params/{name}")
async def params(name: str):
    config = {}
    filepath = f'../../optinist/config/{name}.yaml'
    if os.path.exists(filepath):
        with open(filepath) as f:
            config = yaml.safe_load(f)
    print(config)
    return config
