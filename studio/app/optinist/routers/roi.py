from fastapi import APIRouter, Depends

from studio.app.common.core.workspace.workspace_dependencies import is_workspace_owner
from studio.app.optinist.core.edit_ROI import EditROI, EditRoiUtils
from studio.app.optinist.schemas.roi import RoiList, RoiPos, RoiStatus

router = APIRouter(prefix="/outputs", tags=["outputs"])


@router.post(
    "/image/{filepath:path}/status",
    response_model=RoiStatus,
    dependencies=[Depends(is_workspace_owner)],
)
async def status(filepath: str):
    return EditROI(file_path=filepath).get_status()


@router.post(
    "/image/{filepath:path}/add_roi",
    response_model=bool,
    dependencies=[Depends(is_workspace_owner)],
)
async def add_roi(filepath: str, pos: RoiPos):
    EditROI(file_path=filepath).add(pos)
    return True


@router.post(
    "/image/{filepath:path}/merge_roi",
    response_model=bool,
    dependencies=[Depends(is_workspace_owner)],
)
async def merge_roi(filepath: str, roi_list: RoiList):
    EditROI(file_path=filepath).merge(roi_list.ids)
    return True


@router.post(
    "/image/{filepath:path}/delete_roi",
    response_model=bool,
    dependencies=[Depends(is_workspace_owner)],
)
async def delete_roi(filepath: str, roi_list: RoiList):
    EditROI(file_path=filepath).delete(roi_list.ids)
    return True


@router.post(
    "/image/{filepath:path}/commit_edit",
    response_model=bool,
    dependencies=[Depends(is_workspace_owner)],
)
async def commit_edit(filepath: str):
    EditRoiUtils.execute(filepath)
    return True


@router.post(
    "/image/{filepath:path}/cancel_edit",
    response_model=bool,
    dependencies=[Depends(is_workspace_owner)],
)
async def cancel(filepath: str):
    EditROI(file_path=filepath).cancel()
    return True
