from dataclasses import dataclass
from typing import Dict, List, Union

from pydantic import BaseModel

from studio.app.common.core.snakemake.smk import ForceRun
from studio.app.const import FILETYPE


@dataclass
class NodeType:
    # Data Types
    IMAGE: str = "ImageFileNode"
    CSV: str = "CsvFileNode"
    FLUO: str = "FluoFileNode"
    BEHAVIOR: str = "BehaviorFileNode"
    HDF5: str = "HDF5FileNode"
    MATLAB: str = "MatlabFileNode"
    MICROSCOPE: str = "MicroscopeFileNode"

    # Data Type (Includs above DataType Nodes)
    DATA: str = "DataNode"

    # Algo Type
    ALGO: str = "AlgorithmNode"


class NodeTypeUtil:
    @staticmethod
    def check_nodetype(node_type: str) -> str:
        """
        Check NodeType (DATA or ALGO) from detailed node type
        """
        if node_type in [
            NodeType.IMAGE,
            NodeType.CSV,
            NodeType.FLUO,
            NodeType.BEHAVIOR,
            NodeType.HDF5,
            NodeType.MATLAB,
            NodeType.MICROSCOPE,
        ]:
            return NodeType.DATA
        elif node_type == NodeType.ALGO:
            return NodeType.ALGO
        else:
            None

    @staticmethod
    def check_nodetype_from_filetype(file_type: str) -> str:
        """
        Check NodeType (DATA or ALGO) from file type
        """
        if file_type in [
            FILETYPE.IMAGE,
            FILETYPE.CSV,
            FILETYPE.BEHAVIOR,
            FILETYPE.HDF5,
            FILETYPE.MATLAB,
            FILETYPE.MICROSCOPE,
        ]:
            return NodeType.DATA
        else:
            None


@dataclass
class OutputType:
    IMAGE: str = "images"
    TIMESERIES: str = "timeseries"
    HEATMAP: str = "heatmap"
    ROI: str = "roi"
    SCATTER: str = "scatter"
    BAR: str = "bar"
    HTML: str = "html"
    LINE: str = "line"
    POLAR: str = "polar"
    HISTOGRAM: str = "histogram"
    PIE: str = "pie"


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
    matPath: str = None


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
