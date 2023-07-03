from fastapi import APIRouter

from studio.app.common.core.utils.config_handler import ConfigReader
from studio.app.common.core.utils.filepath_finder import find_param_filepath
from studio.app.optinist.schemas.nwb import NWBParams

router = APIRouter()


@router.get("/nwb", response_model=NWBParams, tags=["params"])
async def get_nwb_params():
    filepath = find_param_filepath("nwb")
    return ConfigReader.read(filepath)
