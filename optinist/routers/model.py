from dataclasses import dataclass
from typing import Dict, List

from pydantic import BaseModel, Field



@dataclass
class Arg:
    name: str
    type: str
    isNone: bool

@dataclass
class Return:
    name: str
    type: str

@dataclass
class Algo:
    args: List[Arg]
    returns: List[Return]
    parameter: str = None
    path: str = None

@dataclass
class TreeNode:
    path: str
    name: str
    isdir: bool
    nodes: List["TreeNode"]

@dataclass
class FILETYPE:
    IMAGE: str = "image"
    CSV: str = "csv"
    HDF5: str = "hdf5"
    BEHAVIOR: str = "behavior"

class DeleteItem(BaseModel):
    uidList: list

class RoiPos(BaseModel):
    posx : int
    posy : int
    sizex : int
    sizey : int
    
class RoiList(BaseModel):
    ids: List[int] = Field(default=[0, 1])
    
class EditRoiSuccess(BaseModel):
    data: List[List] = Field(default=[
        [1,1,1, None, None],
        [None,None,1, None, None],
        [3,3,3, None, None],
        [1,1,1, None, 3],
        [2,2,None, None, 2],
    ])
    max_index: int
@dataclass
class HDF5Node:
    isDir: bool
    name: str
    path: str
    nodes: List['HDF5Node'] = None
    shape: tuple = None
    nbytes: str = None

@dataclass
class OutputData:
    data: Dict[str, dict]
    columns: List[str] = None
    index: List[str] = None

@dataclass
class JsonTimeSeriesData(OutputData):
    xrange: list = None
    std: Dict[str, dict] = None
