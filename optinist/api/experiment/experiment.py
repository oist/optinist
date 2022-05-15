from dataclasses import dataclass
from typing import Dict

from optinist.api.workflow.workflow import Edge, Node


@dataclass
class ExptFunction:
    unique_id: str
    name: str
    success: str


@dataclass
class ExptConfig:
    timestamp: str
    name: str
    unique_id: str
    function: Dict[str, ExptFunction]
    nodeDict: Dict[str, Node]
    edgeDict: Dict[str, Edge]
