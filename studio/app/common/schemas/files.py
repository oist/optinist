from dataclasses import dataclass
from typing import List, Optional

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


@dataclass
class DownloadStatus:
    total: int = 0
    current: int = 0
    error: Optional[str] = None
