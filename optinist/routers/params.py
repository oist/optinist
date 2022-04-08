from fastapi import APIRouter
from optinist.cui_api.config_reader import ConfigReader

router = APIRouter()


@router.get("/params/{name}")
async def params(name: str):
    config = ConfigReader.read(name)
    return config
