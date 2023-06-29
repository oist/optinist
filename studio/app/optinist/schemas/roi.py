from typing import List

from pydantic import BaseModel, Field


class RoiPos(BaseModel):
    posx: int
    posy: int
    sizex: int
    sizey: int


class RoiList(BaseModel):
    ids: List[int] = Field(default=[0, 1])


class EditRoiSuccess(BaseModel):
    max_index: int
