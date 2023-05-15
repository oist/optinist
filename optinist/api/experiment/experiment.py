from dataclasses import dataclass
from typing import Dict, Optional

from optinist.api.workflow.workflow import Edge, Node, OutputPath


@dataclass
class ExptFunction:
    unique_id: str
    name: str
    success: str
    hasNWB: bool
    message: Optional[str] = None
    outputPaths: Optional[Dict[str, OutputPath]] = None
    started_at: Optional[str] = None
    finished_at: Optional[str] = None


@dataclass
class ExptConfig:
    created_at: str
    finished_at: Optional[str]
    name: Optional[str]
    unique_id: str
    hasNWB: bool
    function: Dict[str, ExptFunction]
    nodeDict: Dict[str, Node]
    edgeDict: Dict[str, Edge]


@dataclass
class ExptImportData:
    nodeDict: Dict[str, Node]
    edgeDict: Dict[str, Edge]
