from typing import List, Optional

from pydantic import BaseModel, EmailStr, Field

password_regex = r"^(?=.*\d)(?=.*[!#$%&()*+,-./@_|])(?=.*[a-zA-Z]).{6,255}$"


class User(BaseModel):
    uid: str
    email: EmailStr

    @property
    def is_admin(self) -> bool:
        return False


class ListUserPaging(BaseModel):
    data: Optional[List[User]]
    total_page: Optional[int]


class UserCreate(BaseModel):
    email: EmailStr
    password: str = Field(max_length=255, regex=password_regex)

    class Config:
        anystr_strip_whitespace = True


class UserUpdate(BaseModel):
    email: Optional[EmailStr]


class UserPasswordUpdate(BaseModel):
    old_password: str
    new_password: str = Field(max_length=255, regex=password_regex)

    class Config:
        anystr_strip_whitespace = True
