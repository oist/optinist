from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel

from studio.app.common.schemas.users import UserInfo


class Workspace(BaseModel):
    id: Optional[int]
    name: str
    user: Optional[UserInfo]
    shared_count: int
    created_at: Optional[datetime]
    updated_at: Optional[datetime]

    class Config:
        orm_mode = True


class WorkspaceCreate(BaseModel):
    name: str


class WorkspaceUpdate(WorkspaceCreate):
    pass


class WorkspaceSharePostStatus(BaseModel):
    user_ids: Optional[List[int]]


class WorkspaceShareStatus(BaseModel):
    users: Optional[List[UserInfo]]
