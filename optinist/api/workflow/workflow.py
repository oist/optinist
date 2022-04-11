from dataclasses import dataclass
from typing import Dict


@dataclass
class NodeType:
    IMAGE: str = "ImageFileNode"
    CSV: str = "CsvFileNode"
    FLUO: str = "FluoFileNode"
    BEHAVIOR: str = "BehaviorFileNode"
    HDF5: str = "HDF5FileNode"
    ALGO: str = "AlgorithmNode"


@dataclass
class OutputPath:
    path: str
    type: str
    max_index: int = None


@dataclass
class Message:
    status: str
    message: str
    outputPaths: Dict[str, OutputPath] = None


@dataclass
class ExpFunction:
    unique_id: str
    name: str
    success: str

@dataclass
class NodeData:
    label: str
    param: dict
    path: list or str
    type: str
    fileType: str = None
    hdf5Path: str = None

@dataclass
class NodePosition:
    x: int
    y: int

@dataclass
class Style:
    border: str = None
    height: int = None
    padding: int = None
    width: int = None
    borderRadius: int = None

@dataclass
class Node:
    id: str
    type: str
    data: NodeData
    position: NodePosition
    style: Style

@dataclass
class Edge:
    id: str
    type: str
    animated: bool
    source: str
    sourceHandle: str
    target: str
    targetHandle: str
    style: Style