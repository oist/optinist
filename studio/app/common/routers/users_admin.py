from fastapi import APIRouter, Depends

from studio.app.common.core.auth.auth_dependencies import get_admin_user
from studio.app.common.core.users import crud_users
from studio.app.common.schemas.users import ListUserPaging, User, UserCreate, UserUpdate

router = APIRouter(
    prefix="/admin/users", tags=["admin"], dependencies=[Depends(get_admin_user)]
)


@router.get("/", response_model=ListUserPaging)
async def list_user(offset: int = 0, limit: int = 10):
    return await crud_users.list_user(offset, limit)


@router.post("/", response_model=User)
async def create_user(data: UserCreate):
    return await crud_users.create_user(data)


@router.get("/{user_id}", response_model=User)
async def get_user(user_id: str):
    return await crud_users.get_user(user_id)


@router.put("/{user_id}", response_model=User)
async def update_user(user_id: str, data: UserUpdate):
    return await crud_users.update_user(user_id, data)


@router.delete("/{user_id}")
async def delete_user(user_id: str):
    return await crud_users.delete_user(user_id)
