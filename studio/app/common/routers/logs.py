from typing import List

from fastapi import APIRouter, Depends, HTTPException, Query
from typing_extensions import Optional

from studio.app.common.core.auth.auth_dependencies import get_current_user
from studio.app.common.core.logger import AppLogger
from studio.app.common.core.utils.log_reader import LogLevel, LogReader
from studio.app.common.schemas.outputs import PaginatedLineResult
from studio.app.common.schemas.users import User

router = APIRouter(prefix="/logs", tags=["logs"])

logger = AppLogger.get_logger()


@router.get(
    "",
    summary="Fetch log data with pagination",
)
async def get_log_data(
    current_user: User = Depends(get_current_user),
    offset: int = Query(
        default=-1,
        ge=-1,
        description="The starting position in the log file from which to fetch data."
        "A value of `-1` indicates the request should start from the end of the file",
    ),
    limit: int = Query(
        default=50,
        ge=0,
        description="Max number of log unit to return.",
    ),
    reverse: bool = Query(
        default=True,
        description="Fetch logs in reverse order.",
    ),
    search: Optional[str] = Query(default=None),
    levels: List[LogLevel] = Query(default=[LogLevel.ALL]),
):
    try:
        stop_offset = None
        log_reader = LogReader(
            levels=levels, filter_user_id=current_user.uid if current_user else None
        )

        if search:
            stop_offset, offset = log_reader.get_unit_position_from_search_text(
                search, offset, reverse, search_match_case=False
            )
            if stop_offset is None:
                return PaginatedLineResult(
                    next_offset=offset,
                    prev_offset=offset,
                    data=[],
                )

            logs = log_reader.read_from_offset(
                offset=stop_offset if reverse else offset,
                stop_offset=offset if reverse else stop_offset,
                reverse=False,
            )
            extra_logs = log_reader.read_from_offset(
                offset=logs.prev_offset if reverse else logs.next_offset,
                stop_offset=None,
                limit=limit,
                reverse=reverse,
            )
            return PaginatedLineResult(
                next_offset=max(logs.next_offset, extra_logs.next_offset),
                prev_offset=min(logs.prev_offset, extra_logs.prev_offset),
                data=(
                    extra_logs.data + logs.data
                    if reverse
                    else logs.data + extra_logs.data
                ),
            )

        return log_reader.read_from_offset(
            offset=offset,
            stop_offset=stop_offset,
            limit=limit,
            reverse=reverse,
        )
    except Exception as e:
        logger.error(e, exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))
