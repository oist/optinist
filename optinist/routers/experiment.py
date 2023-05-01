from typing import Dict
from fastapi import APIRouter, HTTPException, status
from fastapi.responses import FileResponse

import shutil
from glob import glob

from optinist.api.config.config_reader import ConfigReader
from optinist.api.config.config_writer import ConfigWriter
from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import join_filepath
from optinist.api.experiment.experiment_reader import ExptConfigReader
from optinist.api.experiment.experiment import ExptConfig, ExptImportData
from optinist.routers.model import DeleteItem, RenameItem


router = APIRouter()


@router.get("/experiments", response_model=Dict[str, ExptConfig], tags=['experiments'])
async def get_experiments():
    exp_config = {}
    config_paths = glob(join_filepath([DIRPATH.OUTPUT_DIR, "*", DIRPATH.EXPERIMENT_YML]))
    for path in config_paths:
        try:
            config = ExptConfigReader.read(path)
            config.nodeDict = []
            config.edgeDict = []
            exp_config[config.unique_id] = config
        except Exception as e:
            pass

    return exp_config


@router.patch("/experiments/{unique_id}/rename", response_model=ExptConfig, tags=['experiments'])
async def rename_experiment(unique_id: str, item: RenameItem):
    config = ExptConfigReader.rename(join_filepath([DIRPATH.OUTPUT_DIR, unique_id, DIRPATH.EXPERIMENT_YML]), new_name=item.new_name)
    config.nodeDict = []
    config.edgeDict = []

    return config


@router.get("/experiments/import/{unique_id}", response_model=ExptImportData, tags=['experiments'])
async def import_experiment(unique_id: str):
    config = ExptConfigReader.read(join_filepath([
        DIRPATH.OUTPUT_DIR,
        unique_id,
        DIRPATH.EXPERIMENT_YML
    ]))
    return {
        "nodeDict": config.nodeDict,
        "edgeDict": config.edgeDict,
    }


def get_last_workflow_uid():
    last_workflow = ConfigReader.read(join_filepath([DIRPATH.OUTPUT_DIR, DIRPATH.WORKFLOW_YML]))
    last_unique_id = last_workflow.get("uid")
    return last_unique_id


def clear_last_workflow_uid():
    ConfigWriter.write(
        dirname=DIRPATH.OUTPUT_DIR,
        filename=DIRPATH.WORKFLOW_YML,
        config={'uid': None}
    )


@router.get("/experiments/last", tags=['experiments'])
async def get_last_experiment() -> str:
    last_unique_id = get_last_workflow_uid()
    if last_unique_id:
        return last_unique_id
    else:
        raise HTTPException(status.HTTP_422_UNPROCESSABLE_ENTITY)


@router.delete("/experiments/{unique_id}", response_model=bool, tags=['experiments'])
async def delete_experiment(unique_id: str):
    try:
        if get_last_workflow_uid() == unique_id:
            clear_last_workflow_uid()
        shutil.rmtree(join_filepath([DIRPATH.OUTPUT_DIR, unique_id]))
        return True
    except Exception as e:
        return False


@router.post("/experiments/delete", response_model=bool, tags=['experiments'])
async def delete_experiment_list(deleteItem: DeleteItem):
    try:
        [
            shutil.rmtree(join_filepath([DIRPATH.OUTPUT_DIR, uid]))
            for uid in deleteItem.uidList
        ]
        if get_last_workflow_uid() in deleteItem.uidList:
            clear_last_workflow_uid()
        return True
    except Exception as e:
        return False


@router.get("/experiments/download/nwb/{unique_id}", tags=['experiments'])
async def download_nwb_experiment(unique_id: str):
    nwb_path_list = glob(join_filepath([
        DIRPATH.OUTPUT_DIR,
        unique_id,
        "*.nwb"
    ]))
    if len(nwb_path_list) > 0:
        return FileResponse(nwb_path_list[0])
    else:
        return False


@router.get("/experiments/download/nwb/{unique_id}/{function_id}", tags=['experiments'])
async def download_nwb_experiment(unique_id: str, function_id: str):
    nwb_path_list = glob(join_filepath([
        DIRPATH.OUTPUT_DIR,
        unique_id,
        function_id,
        "*.nwb"
    ]))
    if len(nwb_path_list) > 0:
        return FileResponse(nwb_path_list[0])
    else:
        return False


@router.get("/experiments/download/config/{unique_id}", tags=['experiments'])
async def download_config_experiment(unique_id: str):
    config_filepath = join_filepath([
        DIRPATH.OUTPUT_DIR,
        unique_id,
        DIRPATH.SNAKEMAKE_CONFIG_YML
    ])
    return FileResponse(config_filepath)
