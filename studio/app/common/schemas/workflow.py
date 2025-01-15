from dataclasses import dataclass
from typing import Dict, Optional

from psutil import Process

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


@dataclass
class WorkflowPIDFileData:
    last_pid: int
    last_script_file: str
    create_time: float
    elapsed_time: str = None

    def __post_init__(self):
        import time

        elapsed_time = time.time() - self.create_time
        self.elapsed_time = f"{elapsed_time:.6}"  # Converted to str for saving to json

    @property
    def elapsed_time_float(self):
        return float(self.elapsed_time)


@dataclass
class WorkflowProcessInfo:
    process: Process
    pid_data: WorkflowPIDFileData


@dataclass
class WorkflowErrorInfo:
    has_error: bool
    error_log: str
