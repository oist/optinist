
from fastapi import APIRouter
from cui_api.get_config_params import get_config_params

router = APIRouter()


@router.get("/nwb")
async def params():
    config = get_config_params("nwb")
    return config
