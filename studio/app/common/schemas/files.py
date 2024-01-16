from dataclasses import dataclass
from typing import List, Optional

from pydantic.dataclasses import dataclass as pydantic_dataclass
from pydantic.networks import AnyHttpUrl


@pydantic_dataclass
class TreeNode:
    path: str
    name: str
    isdir: bool
    nodes: List["TreeNode"]
    shape: Optional[List] = None


@dataclass
class FilePath:
    file_path: str


@dataclass
class DownloadFileRequest:
    url: AnyHttpUrl


@dataclass
class DownloadStatus:
    total: int = 0
    current: int = 0
    error: Optional[str] = None
