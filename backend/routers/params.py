from fastapi import APIRouter
from cui_api.get_config_params import get_config_params

router = APIRouter()


@router.get("/params/{name}")
async def params(name: str):
    config = get_config_params(name)
    return config
