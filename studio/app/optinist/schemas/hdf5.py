from typing import List

from pydantic.dataclasses import dataclass as pydantic_dataclass


@pydantic_dataclass
class HDF5Node:
    isDir: bool
    name: str
    path: str
    nodes: List["HDF5Node"] = None
    shape: tuple = None
    nbytes: str = None
    dataType: str = None
