import shutil
from glob import glob
from typing import Dict

from fastapi import APIRouter
from fastapi.responses import FileResponse

from studio.app.common.core.experiment.experiment import ExptConfig, ExptImportData
from studio.app.common.core.experiment.experiment_reader import ExptConfigReader
from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.schemas.experiment import DeleteItem, RenameItem
from studio.app.dir_path import DIRPATH

router = APIRouter(prefix="/experiments", tags=["experiments"])


@router.get("/{workspace_id}", response_model=Dict[str, ExptConfig])
async def get_experiments(workspace_id: str):
    exp_config = {}
    config_paths = glob(
        join_filepath([DIRPATH.OUTPUT_DIR, workspace_id, "*", DIRPATH.EXPERIMENT_YML])
    )
    for path in config_paths:
        try:
            config = ExptConfigReader.read(path)
            config.nodeDict = []
            config.edgeDict = []
            exp_config[config.unique_id] = config
        except Exception:
            pass

    return exp_config


@router.patch("/{workspace_id}/{unique_id}/rename", response_model=ExptConfig)
async def rename_experiment(workspace_id: str, unique_id: str, item: RenameItem):
    config = ExptConfigReader.rename(
        join_filepath(
            [DIRPATH.OUTPUT_DIR, workspace_id, unique_id, DIRPATH.EXPERIMENT_YML]
        ),
        new_name=item.new_name,
    )
    config.nodeDict = []
    config.edgeDict = []

    return config


@router.get("/import/{workspace_id}/{unique_id}", response_model=ExptImportData)
async def import_experiment(workspace_id: str, unique_id: str):
    config = ExptConfigReader.read(
        join_filepath(
            [DIRPATH.OUTPUT_DIR, workspace_id, unique_id, DIRPATH.EXPERIMENT_YML]
        )
    )
    return {
        "nodeDict": config.nodeDict,
        "edgeDict": config.edgeDict,
    }


@router.delete("/{workspace_id}/{unique_id}", response_model=bool)
async def delete_experiment(workspace_id: str, unique_id: str):
    try:
        shutil.rmtree(join_filepath([DIRPATH.OUTPUT_DIR, workspace_id, unique_id]))
        return True
    except Exception:
        return False


@router.post("/delete/{workspace_id}", response_model=bool)
async def delete_experiment_list(workspace_id: str, deleteItem: DeleteItem):
    try:
        [
            shutil.rmtree(join_filepath([DIRPATH.OUTPUT_DIR, workspace_id, uid]))
            for uid in deleteItem.uidList
        ]
        return True
    except Exception:
        return False


@router.get("/download/config/{workspace_id}/{unique_id}")
async def download_config_experiment(workspace_id: str, unique_id: str):
    config_filepath = join_filepath(
        [DIRPATH.OUTPUT_DIR, workspace_id, unique_id, DIRPATH.SNAKEMAKE_CONFIG_YML]
    )
    return FileResponse(config_filepath)
