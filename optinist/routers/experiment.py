from fastapi import APIRouter
from fastapi.responses import FileResponse

import shutil
from glob import glob
from pydantic import BaseModel

from optinist.cui_api.dir_path import DIRPATH
from optinist.cui_api.filepath_creater import join_filepath
from optinist.cui_api.experiment_config import ExpConfigReader

router = APIRouter()


class DeleteItem(BaseModel):
    uidList: list


@router.get("/experiments")
async def read_experiment():
    exp_config = {}
    config_paths = glob(join_filepath([DIRPATH.BASE_DIR, "*", DIRPATH.EXPERIMENT_YML]))
    for path in config_paths:
        config = ExpConfigReader.read(path)
        config.nodeList = []
        config.edgeList = []
        exp_config[config.unique_id] = config

    return exp_config


@router.get("/experiments/import/{unique_id}")
async def read_experiment(unique_id: str):
    config = ExpConfigReader.read(join_filepath([
        DIRPATH.BASE_DIR, unique_id, DIRPATH.EXPERIMENT_YML]))
    return {
        "nodeList": config.nodeList,
        "edgeList": config.edgeList,
    }


@router.delete("/experiments/{unique_id}")
async def delete_experiment(unique_id: str):
    try:
        shutil.rmtree(join_filepath([DIRPATH.BASE_DIR, unique_id]))
        return True
    except Exception as e:
        return False


@router.post("/experiments/delete")
async def delete_experiment_list(deleteItem: DeleteItem):
    try:
        [shutil.rmtree(join_filepath([DIRPATH.BASE_DIR, uid])) for uid in deleteItem.uidList]
        return True
    except Exception as e:
        return False


@router.get("/experiments/download/{unique_id}")
async def download_experiment(unique_id: str):
    nwb_file = glob(join_filepath([DIRPATH.BASE_DIR, unique_id, "*", "*.nwb"]))[0]
    return FileResponse(nwb_file)