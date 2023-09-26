from typing import Union

from fastapi import Depends, HTTPException
from sqlmodel import Session, or_

from studio.app.common import models as common_model
from studio.app.common.core.auth.auth_dependencies import get_current_user
from studio.app.common.db.database import get_db
from studio.app.common.schemas.users import User


async def is_workspace_owner(
    workspace_id: Union[int, str],
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
):
    workspace = (
        db.query(common_model.Workspace)
        .filter(
            common_model.Workspace.id == workspace_id,
            common_model.Workspace.deleted.is_(False),
            common_model.Workspace.user_id == current_user.id,
        )
        .first()
    )

    if workspace is None:
        raise HTTPException(status_code=403, detail="Operation is not available")
    else:
        return True


async def is_workspace_available(
    workspace_id: Union[int, str],
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> bool:
    workspace = (
        db.query(common_model.Workspace)
        .join(
            common_model.WorkspacesShareUser,
            common_model.WorkspacesShareUser.workspace_id == common_model.Workspace.id,
            isouter=True,
        )
        .filter(
            common_model.Workspace.id == workspace_id,
            common_model.Workspace.deleted.is_(False),
            or_(
                common_model.Workspace.user_id == current_user.id,
                common_model.WorkspacesShareUser.user_id == current_user.id,
            ),
        )
        .group_by(common_model.Workspace.id)
        .first()
    )

    if workspace is None:
        raise HTTPException(status_code=403, detail="Workspace is not available")
    else:
        return True
