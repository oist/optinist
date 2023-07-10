from typing import Optional

from fastapi import Depends, HTTPException, Response, status
from fastapi.security import APIKeyHeader, HTTPAuthorizationCredentials, HTTPBearer
from firebase_admin import auth as firebase_auth

from studio.app.common.core.auth.auth_config import AUTH_CONFIG
from studio.app.common.core.auth.security import validate_access_token
from studio.app.common.core.users.crud_users import get_user
from studio.app.common.schemas.users import User


async def get_current_user(
    res: Response,
    ex_token: Optional[str] = Depends(APIKeyHeader(name="ExToken", auto_error=False)),
    credential: HTTPAuthorizationCredentials = Depends(HTTPBearer(auto_error=False)),
):
    if AUTH_CONFIG.USE_FIREBASE_TOKEN:
        if credential is None:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Unauthorized",
                headers={"WWW-Authenticate": "Bearer realm='auth_required'"},
            )
        try:
            user = firebase_auth.verify_id_token(credential.credentials)
            authed_user: User = await get_user(user["uid"])
            return authed_user
        except Exception:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Could not validate credentials",
                headers={"WWW-Authenticate": "Bearer realm='invalid_token'"},
            )

    if ex_token is None:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            headers={"WWW-Authenticate": 'Bearer realm="auth_required"'},
        )

    payload, err = validate_access_token(ex_token)

    if err:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            headers={"WWW-Authenticate": 'Bearer realm="auth_required"'},
            detail=str(err),
        )

    try:
        return await get_user(payload["sub"])
    except Exception:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            headers={"WWW-Authenticate": 'Bearer realm="auth_required"'},
            detail="Could not validate credentials",
        )


async def get_admin_user(current_user: User = Depends(get_current_user)):
    if current_user.is_admin:
        return current_user
    else:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Insufficient privileges",
        )
