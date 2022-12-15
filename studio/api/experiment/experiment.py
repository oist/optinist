from dataclasses import dataclass
from typing import Dict

from studio.api.workflow.workflow import Edge, Node


@dataclass
class ExptFunction:
    unique_id: str
    name: str
    success: str
    hasNWB: bool


@dataclass
class ExptConfig:
    timestamp: str
    name: str
    unique_id: str
    hasNWB: bool
    function: Dict[str, ExptFunction]
    nodeDict: Dict[str, Node]
    edgeDict: Dict[str, Edge]
