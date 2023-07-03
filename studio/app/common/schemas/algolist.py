from dataclasses import dataclass
from typing import Dict, List, Union

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
