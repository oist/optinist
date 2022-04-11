from dataclasses import dataclass
from typing import Dict, List

from optinist.api.workflow.workflow import Edge, ExpFunction, Node

@dataclass
class ExpConfig:
    timestamp: str
    name: str
    unique_id: str
    function: Dict[str, ExpFunction]
    nodeList: List[Node]
    edgeList: List[Edge]
