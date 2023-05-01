from dataclasses import dataclass
from typing import Dict, List, Union

from pydantic import BaseModel

from optinist.api.snakemake.smk import ForceRun


@dataclass
class NodeType:
    IMAGE: str = "ImageFileNode"
    CSV: str = "CsvFileNode"
    FLUO: str = "FluoFileNode"
    BEHAVIOR: str = "BehaviorFileNode"
    HDF5: str = "HDF5FileNode"
    ALGO: str = "AlgorithmNode"


@dataclass
class OutputType:
    IMAGE: str = "images"
    TIMESERIES: str = "timeseries"
    HEATMAP: str = "heatmap"
    ROI: str = "roi"
    SCATTER: str = "scatter"
    BAR: str = "bar"
    HTML: str = "html"


class NodeItem(BaseModel):
    pendingNodeIdList: list = []


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
class NodeData:
    label: str
    param: dict
    path: Union[str, List]
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


class RunItem(BaseModel):
    name: str = None
    nodeDict: dict = {}
    edgeDict: dict = {}
    snakemakeParam: dict = {}
    nwbParam: dict = {}
    forceRunList: List[ForceRun]
