import os
import shutil
from dataclasses import asdict

import yaml
from fastapi import APIRouter, Depends, File, HTTPException, UploadFile
from fastapi.responses import FileResponse

from studio.app.common.core.experiment.experiment_reader import ExptConfigReader
from studio.app.common.core.experiment.experiment_utils import ExptUtils
from studio.app.common.core.utils.filepath_creater import (
    create_directory,
    join_filepath,
)
from studio.app.common.core.workflow.workflow_reader import WorkflowConfigReader
from studio.app.common.core.workspace.workspace_dependencies import (
    is_workspace_available,
)
from studio.app.common.schemas.workflow import WorkflowConfig, WorkflowWithResults
from studio.app.dir_path import DIRPATH

router = APIRouter(prefix="/workflow", tags=["workflow"])


@router.get(
    "/fetch/{workspace_id}",
    response_model=WorkflowWithResults,
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
        return WorkflowWithResults(
            **asdict(last_expt_config), **asdict(workflow_config)
        )
    else:
        raise HTTPException(status_code=404)


@router.get(
    "/reproduce/{workspace_id}/{unique_id}",
    response_model=WorkflowWithResults,
    dependencies=[Depends(is_workspace_available)],
)
async def reproduce_experiment(workspace_id: str, unique_id: str):
    experiment_config_path = join_filepath(
        [DIRPATH.OUTPUT_DIR, workspace_id, unique_id, DIRPATH.EXPERIMENT_YML]
    )
    workflow_config_path = join_filepath(
        [DIRPATH.OUTPUT_DIR, workspace_id, unique_id, DIRPATH.WORKFLOW_YML]
    )
    if os.path.exists(experiment_config_path) and os.path.exists(workflow_config_path):
        experiment_config = ExptConfigReader.read(experiment_config_path)
        workflow_config = WorkflowConfigReader.read(workflow_config_path)
        return WorkflowWithResults(
            **asdict(experiment_config), **asdict(workflow_config)
        )
    else:
        raise HTTPException(status_code=404, detail="file not found")


@router.get(
    "/download/{workspace_id}/{unique_id}",
    dependencies=[Depends(is_workspace_available)],
)
async def download_workspace_config(workspace_id: str, unique_id: str):
    config_filepath = join_filepath(
        [DIRPATH.OUTPUT_DIR, workspace_id, unique_id, DIRPATH.WORKFLOW_YML]
    )
    if os.path.exists(config_filepath):
        return FileResponse(config_filepath, filename=DIRPATH.WORKFLOW_YML)
    else:
        raise HTTPException(status_code=404, detail="file not found")


@router.post("/import")
async def import_workflow_config(file: UploadFile = File(...)):
    try:
        contents = yaml.safe_load(await file.read())
        if contents is None:
            raise HTTPException(status_code=400, detail="Invalid yaml file")
        return WorkflowConfig(**contents)
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Parsing yaml failed: {str(e)}")


@router.get(
    "/sample_data/{workspace_id}",
    dependencies=[Depends(is_workspace_available)],
)
async def copy_sample_data(workspace_id: str):
    folders = ["input", "output"]

    for folder in folders:
        sample_data_dir = join_filepath([DIRPATH.ROOT_DIR, "sample_data", folder])
        user_dir = join_filepath([DIRPATH.DATA_DIR, folder, workspace_id])

        create_directory(user_dir)
        shutil.copytree(sample_data_dir, user_dir, dirs_exist_ok=True)

    return True
