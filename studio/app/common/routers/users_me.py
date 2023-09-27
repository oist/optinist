from fastapi import APIRouter, Depends
from sqlmodel import Session

from studio.app.common.core.auth.auth_dependencies import get_current_user
from studio.app.common.core.users import crud_users
from studio.app.common.db.database import get_db
from studio.app.common.schemas.users import SelfUserUpdate, User, UserPasswordUpdate

router = APIRouter(prefix="/users/me", tags=["users/me"])


@router.get("", response_model=User)
async def me(current_user: User = Depends(get_current_user)):
    return current_user


@router.put("", response_model=User)
async def update_me(
    data: SelfUserUpdate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
):
    return await crud_users.update_user(
        db, current_user.id, data, organization_id=current_user.organization.id
    )


@router.put("/password", response_model=bool)
async def update_password(
    data: UserPasswordUpdate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
):
    return await crud_users.update_password(
        db, current_user.id, data, organization_id=current_user.organization.id
    )


@router.delete("", response_model=bool)
async def delete_me(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
):
    return await crud_users.delete_user(
        db, current_user.id, organization_id=current_user.organization.id
    )
