
from fastapi import APIRouter
from optinist.cui_api.config_reader import ConfigReader

router = APIRouter()


@router.get("/nwb")
async def params():
    return ConfigReader.read("nwb")
