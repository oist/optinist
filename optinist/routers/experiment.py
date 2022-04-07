from fastapi import APIRouter
from glob import glob
import yaml
import shutil
from pydantic import BaseModel
from typing import List
from fastapi.responses import FileResponse

from optinist.cui_api.const import BASE_DIR
from optinist.cui_api.utils import join_file_path

router = APIRouter()


@router.get("/experiments")
async def read_experiment():
    experiment_config = {}
    config_paths = glob(join_file_path([BASE_DIR, "*", "experiment.yaml"]))
    for path in config_paths:
        with open(path, 'r') as f:
            _config = yaml.safe_load(f)
        experiment_config[_config["unique_id"]] = _config

    return experiment_config


@router.get("/experiments/import/{unique_id}")
async def read_experiment(unique_id: str):
    with open(join_file_path([BASE_DIR, unique_id, "experiment.yaml"])) as f:
        config = yaml.safe_load(f)

    response_config = {}
    response_config["nodeList"] = config["nodeList"]
    response_config["edgeList"] = config["edgeList"]

    return response_config


@router.delete("/experiments/{unique_id}")
async def delete_experiment(unique_id: str):
    try:
        shutil.rmtree(join_file_path([BASE_DIR, unique_id]))
        return True
    except Exception as e:
        return False


class DeleteItem(BaseModel):
    uidList: list

@router.post("/experiments/delete")
async def delete_experiment_list(deleteItem: DeleteItem):
    try:
        [shutil.rmtree(join_file_path([BASE_DIR, uid])) for uid in deleteItem.uidList]
        return True
    except Exception as e:
        return False


@router.get("/experiments/download/{unique_id}")
async def download_experiment(unique_id: str):
    # file_path = "/Users/shogoakiyama/Desktop/ann_0.json"
    nwb_file = glob(join_file_path([BASE_DIR, unique_id, "*", "*.nwb"]))[0]
    return FileResponse(nwb_file)
