from fastapi import APIRouter

from studio.app.common.core.auth import auth
from studio.app.common.schemas.auth import AccessToken, RefreshToken, Token, UserAuth

router = APIRouter(prefix="/auth", tags=["auth"])


@router.post("/login", response_model=Token)
async def login(user_data: UserAuth):
    return await auth.authenticate_user(user_data)


@router.post("/refresh", response_model=AccessToken)
async def refresh(refresh_token: RefreshToken):
    return await auth.refresh_current_user_token(refresh_token.refresh_token)


@router.post("/send_reset_password_mail")
async def send_reset_password_mail(email: str):
    return await auth.send_reset_password_mail(email)
