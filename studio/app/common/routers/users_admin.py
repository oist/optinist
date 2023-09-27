from fastapi import APIRouter, Depends
from fastapi_pagination import LimitOffsetPage
from sqlmodel import Session

from studio.app.common.core.auth.auth_dependencies import get_admin_user
from studio.app.common.core.users import crud_users
from studio.app.common.db.database import get_db
from studio.app.common.schemas.base import SortOptions
from studio.app.common.schemas.users import (
    User,
    UserCreate,
    UserSearchOptions,
    UserUpdate,
)

router = APIRouter(prefix="/admin/users", tags=["users/admin"])


@router.get("", response_model=LimitOffsetPage[User])
async def list_user(
    db: Session = Depends(get_db),
    options: UserSearchOptions = Depends(),
    sortOptions: SortOptions = Depends(),
    current_admin: User = Depends(get_admin_user),
):
    return await crud_users.list_user(
        db,
        organization_id=current_admin.organization.id,
        options=options,
        sortOptions=sortOptions,
    )


@router.post("", response_model=User)
async def create_user(
    data: UserCreate,
    db: Session = Depends(get_db),
    current_admin: User = Depends(get_admin_user),
):
    return await crud_users.create_user(
        db, data, organization_id=current_admin.organization.id
    )


@router.get("/{user_id}", response_model=User)
async def get_user(
    user_id: int,
    db: Session = Depends(get_db),
    current_admin: User = Depends(get_admin_user),
):
    return await crud_users.get_user(
        db, user_id, organization_id=current_admin.organization.id
    )


@router.put("/{user_id}", response_model=User)
async def update_user(
    user_id: int,
    data: UserUpdate,
    db: Session = Depends(get_db),
    current_admin: User = Depends(get_admin_user),
):
    return await crud_users.update_user(
        db, user_id, data, organization_id=current_admin.organization.id
    )


@router.delete("/{user_id}", response_model=bool)
async def delete_user(
    user_id: int,
    db: Session = Depends(get_db),
    current_admin: User = Depends(get_admin_user),
):
    return await crud_users.delete_user(
        db, user_id, organization_id=current_admin.organization.id
    )
