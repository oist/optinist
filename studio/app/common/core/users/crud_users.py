from fastapi import HTTPException
from fastapi_pagination.ext.sqlmodel import paginate
from firebase_admin import auth as firebase_auth
from firebase_admin.auth import UserRecord
from sqlalchemy import func
from sqlalchemy.orm import aliased
from sqlmodel import Session, select

from studio.app.common.core.auth.auth import authenticate_user
from studio.app.common.core.logger import AppLogger
from studio.app.common.core.workspace.workspace_services import WorkspaceService
from studio.app.common.models import Role as RoleModel
from studio.app.common.models import User as UserModel
from studio.app.common.models import UserRole as UserRoleModel
from studio.app.common.models.experiment import ExperimentRecord
from studio.app.common.models.workspace import Workspace
from studio.app.common.schemas.auth import UserAuth
from studio.app.common.schemas.base import SortOptions
from studio.app.common.schemas.users import (
    User,
    UserCreate,
    UserPasswordUpdate,
    UserSearchOptions,
    UserUpdate,
)

logger = AppLogger.get_logger()


async def set_role(db: Session, user_id: int, role_id: int, auto_commit=True):
    db.query(UserRoleModel).filter_by(user_id=user_id).delete(synchronize_session=False)
    role_user = UserRoleModel(user_id=user_id, role_id=role_id)
    db.add(role_user)
    db.flush()
    if auto_commit:
        db.commit()


async def get_user(db: Session, user_id: int, organization_id: int) -> User:
    try:
        data = (
            db.query(UserModel, UserRoleModel.role_id)
            .outerjoin(UserRoleModel, UserModel.id == UserRoleModel.user_id)
            .filter(
                UserModel.id == user_id,
                UserModel.active.is_(True),
                UserModel.organization_id == organization_id,
            )
            .first()
        )
        assert data is not None, "User not found"
        user, role_id = data
        user.__dict__["role_id"] = role_id
        return User.from_orm(user)
    except AssertionError as e:
        logger.error(e, exc_info=True)
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(e, exc_info=True)
        raise HTTPException(status_code=400, detail=str(e))


async def list_user(
    db: Session,
    organization_id: int,
    options: UserSearchOptions,
    sortOptions: SortOptions,
):
    def user_transformer(items):
        users = []
        for item in items:
            user, role_id, data_usage = item
            user.__dict__["role_id"] = role_id
            user.__dict__["data_usage"] = data_usage
            users.append(user)
        return users

    try:
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

        WorkspaceCapacity = aliased(workspace_capacity_subq)
        ExperimentCapacity = aliased(experiment_capacity_subq)

        sa_sort_list = sortOptions.get_sa_sort_list(
            sa_table=UserModel,
            mapping={"role_id": UserRoleModel.role_id, "role": RoleModel.role},
        )
        users = paginate(
            db,
            query=select(
                UserModel,
                func.min(UserRoleModel.role_id),
                func.coalesce(WorkspaceCapacity.c.input_workspace_capacity, 0)
                + func.coalesce(ExperimentCapacity.c.experiment_capacity, 0).label(
                    "data_usage"
                ),
            )
            .outerjoin(WorkspaceCapacity, WorkspaceCapacity.c.user_id == UserModel.id)
            .outerjoin(ExperimentCapacity, ExperimentCapacity.c.user_id == UserModel.id)
            .join(UserRoleModel, UserRoleModel.user_id == UserModel.id, isouter=True)
            .join(RoleModel, RoleModel.id == UserRoleModel.role_id, isouter=True)
            .filter(
                UserModel.active.is_(True),
                UserModel.organization_id == organization_id,
            )
            .filter(
                UserModel.name.like("%{0}%".format(options.name)),
                UserModel.email.like("%{0}%".format(options.email)),
            )
            .group_by(UserModel.id)
            .order_by(*sa_sort_list),
            transformer=user_transformer,
            unique=False,
        )
        return users
    except Exception as e:
        logger.error(e, exc_info=True)
        raise HTTPException(status_code=400, detail=str(e))


async def create_user(db: Session, data: UserCreate, organization_id: int):
    try:
        user: UserRecord = firebase_auth.create_user(
            email=data.email, password=data.password
        )
        user_db = UserModel(
            uid=user.uid,
            email=user.email,
            name=data.name,
            organization_id=organization_id,
            active=True,
        )
        # it may be possible to specify the organization_id externally in the future
        db.add(user_db)
        db.flush()
        await set_role(db, user_id=user_db.id, role_id=data.role_id, auto_commit=False)
        db.commit()
        user_db.__dict__["role_id"] = data.role_id
        return User.from_orm(user_db)
    except Exception as e:
        logger.error(e, exc_info=True)
        raise HTTPException(status_code=400, detail=str(e))


async def update_user(
    db: Session, user_id: int, data: UserUpdate, organization_id: int
):
    try:
        user_db = (
            db.query(UserModel)
            .filter(
                UserModel.active.is_(True),
                UserModel.id == user_id,
                UserModel.organization_id == organization_id,
            )
            .first()
        )
        assert user_db is not None, "User not found"
        user_data = data.dict(exclude_unset=True)
        role_id = user_data.pop("role_id", None)
        for key, value in user_data.items():
            setattr(user_db, key, value)
        if role_id is not None:
            await set_role(db, user_id=user_db.id, role_id=role_id, auto_commit=False)
        user_db.__dict__["role_id"] = role_id
        firebase_auth.update_user(user_db.uid, email=data.email)
        db.commit()
        return User.from_orm(user_db)
    except AssertionError as e:
        logger.error(e, exc_info=True)
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(e, exc_info=True)
        raise HTTPException(status_code=400, detail=str(e))


async def update_password(
    db: Session,
    user_id: int,
    data: UserPasswordUpdate,
    organization_id: int,
):
    user = await get_user(db, user_id, organization_id)
    await authenticate_user(
        db, data=UserAuth(email=user.email, password=data.old_password)
    )
    try:
        user = firebase_auth.update_user(user.uid, password=data.new_password)
        return True
    except Exception as e:
        logger.error(e, exc_info=True)
        raise HTTPException(status_code=400, detail=str(e))


async def delete_user(db: Session, user_id: int, organization_id: int) -> bool:
    try:
        # ----------------------------------------
        # Delete a User workspace contents
        # ----------------------------------------

        workspaces = (
            db.query(Workspace)
            .filter(
                Workspace.user_id == user_id,
                Workspace.deleted.is_(False),
            )
            .all()
        )
        workspace_ids = [ws.id for ws in workspaces]

        # Delete owned workspaces
        for workspace_id in workspace_ids:
            WorkspaceService.process_workspace_deletion(
                db=db, workspace_id=workspace_id, user_id=user_id
            )

        # ----------------------------------------
        # Delete a User database record
        # ----------------------------------------

        user_db = (
            db.query(UserModel)
            .filter(
                UserModel.active.is_(True),
                UserModel.id == user_id,
                UserModel.organization_id == organization_id,
            )
            .first()
        )
        assert user_db is not None, "User not found"

        user_db.active = False

        # ----------------------------------------
        # Delete a User firebase account
        # ----------------------------------------

        firebase_auth.delete_user(user_db.uid)

        # The transaction is committed at this point
        # ATTENTION:
        #   - If an exception occurs when deleting a Firebase account,
        #     this commit may not be executed and the account may become undeletable.
        #   - One possible solution to this issue is to add a status
        #     when an error occurs (such as "Account suspended").
        db.commit()

        return True

    except AssertionError as e:
        logger.error(e, exc_info=True)
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(e, exc_info=True)
        raise HTTPException(status_code=400, detail=str(e))
