from typing import List

from pydantic.dataclasses import dataclass as pydantic_dataclass


@pydantic_dataclass
class MatNode:
    isDir: bool
    name: str
    path: str
    nodes: List["MatNode"] = None
    shape: tuple = None
    nbytes: str = None
    dataType: str = None
