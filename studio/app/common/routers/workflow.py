import yaml
from fastapi import APIRouter, File, HTTPException, UploadFile
from fastapi.responses import FileResponse

from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.workflow.workflow import WorkflowConfig
from studio.app.common.core.workflow.workflow_reader import WorkflowConfigReader
from studio.app.dir_path import DIRPATH

router = APIRouter(prefix="/workflow", tags=["workflow"])


@router.get("/import/{workspace_id}/{unique_id}", response_model=WorkflowConfig)
async def import_experiment(workspace_id: str, unique_id: str):
    config = WorkflowConfigReader.read(
        join_filepath(
            [DIRPATH.OUTPUT_DIR, workspace_id, unique_id, DIRPATH.WORKFLOW_YML]
        )
    )
    return config


@router.get("/download/{workspace_id}/{unique_id}")
async def download_workspace_config(workspace_id: str, unique_id: str):
    config_filepath = join_filepath(
        [DIRPATH.OUTPUT_DIR, workspace_id, unique_id, DIRPATH.WORKFLOW_YML]
    )
    return FileResponse(config_filepath, filename=DIRPATH.WORKFLOW_YML)


@router.post("/load")
async def load_workflow_config(file: UploadFile = File(...)):
    try:
        contents = yaml.safe_load(await file.read())
        if contents is None:
            raise HTTPException(status_code=400, detail="Invalid yaml file")
        return WorkflowConfig(**contents)
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Parsing yaml failed: {str(e)}")
