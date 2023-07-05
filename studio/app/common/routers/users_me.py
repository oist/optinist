from fastapi import APIRouter, Depends

from studio.app.common.core.auth.auth_dependencies import get_current_user
from studio.app.common.core.users import crud_users
from studio.app.common.schemas.users import User, UserPasswordUpdate, UserUpdate

router = APIRouter(prefix="/users/me", tags=["mypage"])


@router.get("/", response_model=User)
async def me(current_user: User = Depends(get_current_user)):
    return current_user


@router.put("/", response_model=User)
async def update_me(data: UserUpdate, current_user: User = Depends(get_current_user)):
    return await crud_users.update_user(current_user.uid, data)


@router.put("/password", response_model=User)
async def update_password(
    data: UserPasswordUpdate,
    current_user: User = Depends(get_current_user),
):
    return await crud_users.upate_password(current_user.uid, data)
