from dataclasses import asdict
from typing import Dict

from studio.app.common.core.utils.config_handler import ConfigWriter
from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.workflow.workflow import Edge, Node
from studio.app.common.core.workflow.workflow_builder import WorkflowConfigBuilder
from studio.app.common.schemas.workflow import WorkflowConfig
from studio.app.dir_path import DIRPATH


class WorkflowConfigWriter:
    def __init__(
        self,
        workspace_id: str,
        unique_id: str,
        nodeDict: Dict[str, Node],
        edgeDict: Dict[str, Edge],
    ) -> None:
        self.workspace_id = workspace_id
        self.unique_id = unique_id
        self.nodeDict = nodeDict
        self.edgeDict = edgeDict
        self.builder = WorkflowConfigBuilder()
        self.create_config()

    @staticmethod
    def write_raw(workspace_id: str, unique_id: str, config: dict) -> None:
        ConfigWriter.write(
            dirname=join_filepath([DIRPATH.OUTPUT_DIR, workspace_id, unique_id]),
            filename=DIRPATH.WORKFLOW_YML,
            config=config,
        )

    def write(self) -> None:
        self.write_raw(
            self.workspace_id, self.unique_id, config=asdict(self.builder.build())
        )

    def create_config(self) -> WorkflowConfig:
        return (
            self.builder.set_node_dict(self.nodeDict)
            .set_edge_dict(self.edgeDict)
            .build()
        )
