from math import ceil

from fastapi import HTTPException
from firebase_admin import auth as firebase_auth
from firebase_admin.auth import UserRecord

from studio.app.common.core.auth.auth import authenticate_user
from studio.app.common.schemas.users import (
    ListUserPaging,
    User,
    UserCreate,
    UserPasswordUpdate,
    UserUpdate,
)


async def get_user(uid: str) -> User:
    try:
        user: UserRecord = firebase_auth.get_user(uid)
        return User(uid=user.uid, email=user.email)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


async def list_user(offset: int = 0, limit: int = 10) -> ListUserPaging:
    try:
        total_user = []
        page = firebase_auth.list_users()
        while page:
            [total_user.append(user) for user in page.users]
            page = page.get_next_page()

        users = total_user[
            offset : offset + limit if offset + limit <= len(total_user) else None
        ]
        data = [User(uid=user.uid, email=user.email) for user in users]

        total_page = ceil(len(total_user) / limit)

        return ListUserPaging(data=data, total_page=total_page)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


async def create_user(data: UserCreate):
    try:
        user = firebase_auth.create_user(email=data.email, password=data.password)
        return User(uid=user.uid, email=user.email)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


async def update_user(uid: str, data: UserUpdate):
    try:
        user = firebase_auth.update_user(uid, email=data.email)
        return User(uid=user.uid, email=user.email)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


async def upate_password(uid: str, data: UserPasswordUpdate):
    user = await get_user(uid)
    _, err = await authenticate_user(user.email, data.old_password)

    if err:
        raise HTTPException(status_code=400, detail="Invalid password")
    try:
        user = firebase_auth.update_user(uid, password=data.new_password)
        return User(uid=user.uid, email=user.email)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


async def delete_user(uid: str):
    try:
        user = await get_user(uid)
        firebase_auth.delete_user(user.uid)
        return user
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
