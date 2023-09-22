import os
import shutil
from dataclasses import asdict
from glob import glob
from typing import Dict

from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import FileResponse

from studio.app.common.core.experiment.experiment import ExptConfig
from studio.app.common.core.experiment.experiment_reader import ExptConfigReader
from studio.app.common.core.experiment.experiment_utils import ExptUtils
from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.workflow.workflow_reader import WorkflowConfigReader
from studio.app.common.core.workspace.workspace_dependencies import (
    is_workspace_available,
    is_workspace_owner,
)
from studio.app.common.schemas.experiment import (
    DeleteItem,
    FetchExptResponse,
    RenameItem,
)
from studio.app.dir_path import DIRPATH

router = APIRouter(prefix="/experiments", tags=["experiments"])


@router.get(
    "/{workspace_id}",
    response_model=Dict[str, ExptConfig],
    dependencies=[Depends(is_workspace_available)],
)
async def get_experiments(workspace_id: str):
    exp_config = {}
    config_paths = glob(
        join_filepath([DIRPATH.OUTPUT_DIR, workspace_id, "*", DIRPATH.EXPERIMENT_YML])
    )
    for path in config_paths:
        try:
            config = ExptConfigReader.read(path)
            exp_config[config.unique_id] = config
        except Exception:
            pass

    return exp_config


@router.patch(
    "/{workspace_id}/{unique_id}/rename",
    response_model=ExptConfig,
    dependencies=[Depends(is_workspace_owner)],
)
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


@router.delete(
    "/{workspace_id}/{unique_id}",
    response_model=bool,
    dependencies=[Depends(is_workspace_owner)],
)
async def delete_experiment(workspace_id: str, unique_id: str):
    try:
        shutil.rmtree(join_filepath([DIRPATH.OUTPUT_DIR, workspace_id, unique_id]))
        return True
    except Exception:
        return False


@router.post(
    "/delete/{workspace_id}",
    response_model=bool,
    dependencies=[Depends(is_workspace_owner)],
)
async def delete_experiment_list(workspace_id: str, deleteItem: DeleteItem):
    try:
        [
            shutil.rmtree(join_filepath([DIRPATH.OUTPUT_DIR, workspace_id, uid]))
            for uid in deleteItem.uidList
        ]
        return True
    except Exception:
        return False


@router.get(
    "/fetch/{workspace_id}",
    response_model=FetchExptResponse,
    dependencies=[Depends(is_workspace_available)],
)
async def fetch_last_experiment(workspace_id: str):
    last_expt_config = ExptUtils.get_last_experiment(workspace_id)
    if last_expt_config:
        workflow_config_path = join_filepath(
            [
                DIRPATH.OUTPUT_DIR,
                workspace_id,
                last_expt_config.unique_id,
                DIRPATH.WORKFLOW_YML,
            ]
        )
        workflow_config = WorkflowConfigReader.read(workflow_config_path)
        return FetchExptResponse(**asdict(last_expt_config), **asdict(workflow_config))
    else:
        raise HTTPException(status_code=404)


@router.get(
    "/download/config/{workspace_id}/{unique_id}",
    dependencies=[Depends(is_workspace_available)],
)
async def download_config_experiment(workspace_id: str, unique_id: str):
    config_filepath = join_filepath(
        [DIRPATH.OUTPUT_DIR, workspace_id, unique_id, DIRPATH.SNAKEMAKE_CONFIG_YML]
    )
    if os.path.exists(config_filepath):
        return FileResponse(config_filepath)
    else:
        raise HTTPException(status_code=404, detail="file not found")
