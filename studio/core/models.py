from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Union

from pydantic import BaseModel, Field
from pydantic.dataclasses import dataclass as pydantic_dataclass


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


@pydantic_dataclass
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
    posx: int
    posy: int
    sizex: int
    sizey: int


class RoiList(BaseModel):
    ids: List[int] = Field(default=[0, 1])


class EditRoiSuccess(BaseModel):
    max_index: int


@pydantic_dataclass
class HDF5Node:
    isDir: bool
    name: str
    path: str
    nodes: List["HDF5Node"] = None
    shape: tuple = None
    nbytes: str = None


@dataclass
class OutputData:
    data: Union[List, Dict, str]
    columns: List[str] = None
    index: List[str] = None


@dataclass
class JsonTimeSeriesData(OutputData):
    xrange: list = None
    std: Dict[str, dict] = None


@dataclass
class FilePath:
    file_path: str


class AlgoModel(BaseModel):
    children: Union[Dict[str, Algo], Dict[str, "AlgoModel"]]


class AlgoList(BaseModel):
    __root__: Dict[str, AlgoModel] = {
        "caiman": {
            "children": {
                "caiman_mc": {
                    "args": [{"name": "image", "type": "ImageData", "isNone": False}],
                    "returns": [{"name": "mc_images", "type": "ImageData"}],
                    "parameter": None,
                    "path": "caiman/caiman_mc",
                }
            }
        }
    }


class NWBParams(BaseModel):
    session_description: str = "optinist"
    identifier: str = "optinist"
    experiment_description: Optional[str] = None
    device: Union[Dict, Any]
    optical_channel: Union[Dict, Any]
    imaging_plane: Union[Dict, Any]
    image_series: Union[Dict, Any]
    ophys: Union[Dict, Any]


class SnakemakeParams(BaseModel):
    use_conda: bool
    cores: int
    forceall: bool
    forcetargets: bool
    lock: bool


class RenameItem(BaseModel):
    new_name: str
