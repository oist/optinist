from dataclasses import dataclass
from typing import List

from pydantic.dataclasses import dataclass as pydantic_dataclass


@pydantic_dataclass
class TreeNode:
    path: str
    name: str
    isdir: bool
    nodes: List["TreeNode"]


@dataclass
class FilePath:
    file_path: str
