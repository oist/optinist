import json
import logging

from fastapi import HTTPException, status
from fastapi.responses import JSONResponse
from requests.exceptions import HTTPError
from sqlmodel import Session

from studio.app.common.core.auth import pyrebase_app
from studio.app.common.core.auth.security import (
    create_access_token,
    create_refresh_token,
    validate_refresh_token,
)
from studio.app.common.models.user import User as UserModel
from studio.app.common.schemas.auth import AccessToken, Token, UserAuth


async def authenticate_user(db: Session, data: UserAuth):
    try:
        user = pyrebase_app.auth().sign_in_with_email_and_password(
            data.email, data.password
        )
        user_db = (
            db.query(UserModel)
            .filter(UserModel.uid == user["localId"], UserModel.active.is_(True))
            .first()
        )
        assert user_db is not None, "Invalid user uid"
        ex_token = create_access_token(subject=user_db.uid)
        return Token(
            access_token=user["idToken"],
            refresh_token=create_refresh_token(subject=user["refreshToken"]),
            token_type="bearer",
            ex_token=ex_token,
        )
    except Exception as e:
        logging.getLogger().error(e)
        raise HTTPException(status_code=400, detail=str(e))


async def refresh_current_user_token(refresh_token: str):
    token, err = validate_refresh_token(refresh_token)

    if err:
        raise HTTPException(status_code=400)
    try:
        user = pyrebase_app.auth().refresh(refresh_token=token["sub"])
        return AccessToken(access_token=user["idToken"])
    except Exception as e:
        logging.getLogger().error(e)
        raise HTTPException(status_code=400)


async def send_reset_password_mail(db: Session, email: str):
    try:
        db.query(UserModel).filter(
            UserModel.email == email, UserModel.active.is_(True)
        ).one()
        pyrebase_app.auth().send_password_reset_email(email)
        return JSONResponse(content=None, status_code=status.HTTP_200_OK)
    except HTTPError as e:
        logging.getLogger().error(e)
        err = json.loads(e.strerror)
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=err.get("error").get("message"),
        )
    except Exception as e:
        logging.getLogger().error(e)
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR)
