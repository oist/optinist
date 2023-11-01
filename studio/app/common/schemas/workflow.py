from dataclasses import dataclass
from typing import Dict, Optional

from studio.app.common.core.experiment.experiment import ExptFunction
from studio.app.common.core.workflow.workflow import Edge, Node


@dataclass
class WorkflowConfig:
    nodeDict: Dict[str, Node]
    edgeDict: Dict[str, Edge]


@dataclass
class WorkflowWithResults:
    workspace_id: str
    unique_id: str
    name: str
    started_at: str
    finished_at: Optional[str]
    success: Optional[str]
    hasNWB: bool
    function: Dict[str, ExptFunction]
    nodeDict: Dict[str, Node]
    edgeDict: Dict[str, Edge]
