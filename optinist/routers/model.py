from dataclasses import dataclass
from typing import Dict, List

from pydantic import BaseModel


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

class RenameItem(BaseModel):
    new_name: str

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
