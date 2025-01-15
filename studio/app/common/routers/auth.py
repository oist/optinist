from fastapi import APIRouter, Depends, HTTPException, status
from sqlmodel import Session

from studio.app.common.core.auth import auth
from studio.app.common.core.logger import AppLogger
from studio.app.common.db.database import get_db
from studio.app.common.schemas.auth import AccessToken, RefreshToken, Token, UserAuth

router = APIRouter(prefix="/auth", tags=["auth"])

logger = AppLogger.get_logger()


@router.post("/login", response_model=Token)
async def login(user_data: UserAuth, db: Session = Depends(get_db)):
    try:
        token, user = await auth.authenticate_user(db, user_data)

    except HTTPException as e:
        logger.error(e, exc_info=True)
        raise e

    except Exception as e:
        logger.error(e, exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Some error occurred during authentication.",
        )

    return token


@router.post("/refresh", response_model=AccessToken)
async def refresh(refresh_token: RefreshToken):
    return await auth.refresh_current_user_token(refresh_token.refresh_token)


@router.post("/send_reset_password_mail")
async def send_reset_password_mail(email: str, db: Session = Depends(get_db)):
    return await auth.send_reset_password_mail(db, email)
