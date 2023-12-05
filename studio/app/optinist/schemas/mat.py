from typing import List, Union

from pydantic.dataclasses import dataclass as pydantic_dataclass

from studio.app.common.schemas.outputs import OutputData


@pydantic_dataclass
class MatNode:
    isDir: bool
    name: str
    path: str
    nodes: List["MatNode"] = None
    shape: tuple = None
    nbytes: str = None
    data: Union[int, str, OutputData] = None
