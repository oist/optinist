from typing import List, Union

from pydantic import BaseModel, Field


class RoiPos(BaseModel):
    posx: int
    posy: int
    sizex: int
    sizey: int


class RoiList(BaseModel):
    ids: List[int] = Field(default=[0, 1])


class RoiStatus(BaseModel):
    temp_add_roi: List[Union[int, float]]
    temp_merge_roi: List[Union[int, float]]
    temp_delete_roi: List[Union[int, float]]
