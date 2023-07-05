from typing import Optional

from pydantic import BaseModel, EmailStr


class Token(BaseModel):
    access_token: str
    token_type: str
    refresh_token: Optional[str]
    ex_token: Optional[str]


class RefreshToken(BaseModel):
    refresh_token: str


class UserAuth(BaseModel):
    email: EmailStr
    password: str


class AccessToken(BaseModel):
    access_token: str
