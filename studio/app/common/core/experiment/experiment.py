from dataclasses import dataclass
from typing import Dict, Optional

from studio.app.common.core.workflow.workflow import Edge, Node


@dataclass
class ExptFunction:
    unique_id: str
    name: str
    success: str
    hasNWB: bool


@dataclass
class ExptConfig:
    timestamp: str
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
