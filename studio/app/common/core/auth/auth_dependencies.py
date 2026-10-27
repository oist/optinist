import logging
from typing import Optional

import sqlalchemy
from fastapi import Depends, HTTPException, Response, status
from fastapi.security import APIKeyHeader, HTTPAuthorizationCredentials, HTTPBearer
from pydantic import ValidationError
from sqlalchemy import func
from sqlalchemy.orm import aliased
from sqlmodel import Session, select

from studio.app.common.core.auth.auth_config import AUTH_CONFIG
from studio.app.common.core.auth.auth_helper import (
    extract_uid_from_firebase_credential,
    extract_uid_from_jwt_token,
)
from studio.app.common.db.database import get_db
from studio.app.common.models import User as UserModel
from studio.app.common.models import UserRole as UserRoleModel
from studio.app.common.models.experiment import ExperimentRecord
from studio.app.common.models.workspace import Workspace
from studio.app.common.schemas.users import User


async def get_current_user(
    res: Response,
    ex_token: Optional[str] = Depends(APIKeyHeader(name="ExToken", auto_error=False)),
    credential: HTTPAuthorizationCredentials = Depends(HTTPBearer(auto_error=False)),
    db: Session = Depends(get_db),
) -> User:
    use_firebase_auth = AUTH_CONFIG.USE_FIREBASE_TOKEN
    try:
        assert credential is not None if use_firebase_auth else True
        assert ex_token is not None if not use_firebase_auth else True

        # Extract uid using helper functions
        uid = None
        err = None
        if use_firebase_auth:
            uid, err = extract_uid_from_firebase_credential(credential)
        else:
            uid, err = extract_uid_from_jwt_token(ex_token)

        assert err is None, str(err)
        assert uid is not None, "Failed to extract user ID"

        # Query user record
        user_data = __get_current_user_record(db, uid)
        assert user_data is not None, "Invalid user data"
        authed_user, role_id, data_usage = user_data
        authed_user.__dict__["role_id"] = role_id
        authed_user.__dict__["data_usage"] = data_usage

        return User.from_orm(authed_user)

    except ValidationError as e:
        logging.getLogger().error(e)
        raise HTTPException(status_code=422, detail=f"Validator Error: {e}")
    except Exception as e:
        logging.getLogger().error(e)
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            headers={"WWW-Authenticate": 'Bearer realm="auth_required"'},
            detail=str(e) or "Could not validate credentials",
        )


def __get_current_user_record(db: Session, uid: str) -> sqlalchemy.engine.row.Row:
    # Make query (calc workspace capacity)
    workspace_capacity_subq = (
        select(
            Workspace.user_id,
            func.coalesce(func.sum(Workspace.input_data_usage), 0).label(
                "input_workspace_capacity"
            ),
        )
        .where(Workspace.deleted.is_(False))
        .group_by(Workspace.user_id)
        .subquery()
    )
    WorkspaceCapacity = aliased(workspace_capacity_subq)

    # Make query (calc experient capacity)
    experiment_capacity_subq = (
        select(
            Workspace.user_id,
            func.coalesce(func.sum(ExperimentRecord.data_usage), 0).label(
                "experiment_capacity"
            ),
        )
        .join(ExperimentRecord, ExperimentRecord.workspace_id == Workspace.id)
        .where(Workspace.deleted.is_(False))
        .group_by(Workspace.user_id)
        .subquery()
    )
    ExperimentCapacity = aliased(experiment_capacity_subq)

    # Query user record
    user_data: sqlalchemy.engine.row.Row = (
        db.query(
            UserModel,
            func.min(UserRoleModel.role_id),
            func.coalesce(WorkspaceCapacity.c.input_workspace_capacity, 0)
            + func.coalesce(ExperimentCapacity.c.experiment_capacity, 0).label(
                "data_usage"
            ),
        )
        .outerjoin(WorkspaceCapacity, WorkspaceCapacity.c.user_id == UserModel.id)
        .outerjoin(ExperimentCapacity, ExperimentCapacity.c.user_id == UserModel.id)
        .outerjoin(UserRoleModel, UserRoleModel.user_id == UserModel.id)
        .filter(UserModel.uid == uid)
        .first()
    )

    return user_data


async def get_admin_user(current_user: User = Depends(get_current_user)) -> User:
    if current_user.is_admin:
        return current_user
    else:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Insufficient privileges",
        )
