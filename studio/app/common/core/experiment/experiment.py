from dataclasses import dataclass
from typing import Dict, Optional

from studio.app.common.core.snakemake.smk import SmkParam
from studio.app.common.core.workflow.workflow import OutputPath
from studio.app.optinist.schemas.nwb import NWBParams


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
    workspace_id: str
    unique_id: str
    name: str
    started_at: str
    finished_at: Optional[str]
    success: Optional[str]
    hasNWB: bool
    function: Dict[str, ExptFunction]
    nwb: NWBParams
    snakemake: SmkParam


@dataclass
class ExptOutputPathIds:
    output_dir: Optional[str] = None
    workspace_id: Optional[str] = None
    unique_id: Optional[str] = None
    function_id: Optional[str] = None

    def __post_init__(self):
        """
        Extract each ID from output_path
        - output_dir format
          - {DIRPATH.OUTPUT_DIR}/{workspace_id}/{unique_id}/{function_id}
        """
        if self.output_dir:
            splitted_path = self.output_dir.split("/")
        else:
            splitted_path = []

        if len(splitted_path) >= 3:
            self.workspace_id, self.unique_id, self.function_id = splitted_path[-3:]
