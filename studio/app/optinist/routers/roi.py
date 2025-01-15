from fastapi import APIRouter, Depends, HTTPException, status

from studio.app.common.core.logger import AppLogger
from studio.app.common.core.workspace.workspace_dependencies import is_workspace_owner
from studio.app.optinist.core.edit_ROI import EditROI, EditRoiUtils
from studio.app.optinist.schemas.roi import RoiList, RoiPos, RoiStatus

router = APIRouter(prefix="/outputs", tags=["outputs"])

logger = AppLogger.get_logger()


@router.post(
    "/image/{filepath:path}/status",
    response_model=RoiStatus,
    dependencies=[Depends(is_workspace_owner)],
)
async def status_roi(filepath: str):
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
    try:
        EditRoiUtils.execute(filepath)

    except Exception as e:
        logger.error(e, exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to commit Edit ROI",
        )

    return True


@router.post(
    "/image/{filepath:path}/cancel_edit",
    response_model=bool,
    dependencies=[Depends(is_workspace_owner)],
)
async def cancel_edit(filepath: str):
    EditROI(file_path=filepath).cancel()
    return True
