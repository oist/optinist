from glob import glob

from fastapi import APIRouter, Depends
from fastapi.responses import FileResponse

from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.workflow.workflow_params import read_default_params
from studio.app.common.core.workspace.workspace_dependencies import (
    is_workspace_available,
)
from studio.app.dir_path import DIRPATH
from studio.app.optinist.schemas.nwb import NWBParams

router = APIRouter()


@router.get("/nwb", response_model=NWBParams, tags=["params"])
async def get_nwb_params():
    return read_default_params("nwb")


@router.get(
    "/experiments/download/nwb/{workspace_id}/{unique_id}",
    dependencies=[Depends(is_workspace_available)],
    tags=["experiments"],
)
async def download_nwb_experiment(workspace_id: str, unique_id: str):
    nwb_path_list = glob(
        join_filepath([DIRPATH.OUTPUT_DIR, workspace_id, unique_id, "*.nwb"])
    )
    if len(nwb_path_list) > 0:
        return FileResponse(nwb_path_list[0])
    else:
        return False


@router.get(
    "/experiments/download/nwb/{workspace_id}/{unique_id}/{function_id}",
    dependencies=[Depends(is_workspace_available)],
    tags=["experiments"],
)
async def download_nwb_experiment_with_function_id(
    workspace_id: str, unique_id: str, function_id: str
):
    nwb_path_list = glob(
        join_filepath(
            [DIRPATH.OUTPUT_DIR, workspace_id, unique_id, function_id, "*.nwb"]
        )
    )
    if len(nwb_path_list) > 0:
        return FileResponse(nwb_path_list[0])
    else:
        return False
