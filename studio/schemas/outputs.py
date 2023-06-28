from dataclasses import dataclass
from typing import Dict, List, Union

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


@dataclass
class OutputData:
    data: Union[List, Dict, str]
    columns: List[str] = None
    index: List[str] = None


@dataclass
class JsonTimeSeriesData(OutputData):
    xrange: list = None
    std: Dict[str, dict] = None
