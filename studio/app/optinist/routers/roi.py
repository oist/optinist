from fastapi import APIRouter, Depends

from studio.app.common.core.workspace.workspace_dependencies import is_workspace_owner
from studio.app.optinist.core.edit_ROI import ACTION, EditROI
from studio.app.optinist.schemas.roi import EditRoiSuccess, RoiList, RoiPos

router = APIRouter(prefix="/outputs", tags=["outputs"])


@router.post(
    "/image/{filepath:path}/add_roi",
    response_model=EditRoiSuccess,
    dependencies=[Depends(is_workspace_owner)],
)
async def add_roi(filepath: str, pos: RoiPos):
    EditROI(action=ACTION.ADD, filepath=filepath, params=pos.dict()).excute()
    return EditRoiSuccess(max_index=0)


@router.post(
    "/image/{filepath:path}/merge_roi",
    response_model=EditRoiSuccess,
    dependencies=[Depends(is_workspace_owner)],
)
async def merge_roi(filepath: str, roi_list: RoiList):
    EditROI(action=ACTION.MERGE, filepath=filepath, params=roi_list.dict()).excute()
    return EditRoiSuccess(max_index=0)


@router.post(
    "/image/{filepath:path}/delete_roi",
    response_model=EditRoiSuccess,
    dependencies=[Depends(is_workspace_owner)],
)
async def delete_roi(filepath: str, roi_list: RoiList):
    EditROI(action=ACTION.DELETE, filepath=filepath, params=roi_list.dict()).excute()
    return EditRoiSuccess(max_index=0)
