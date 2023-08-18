import os
from dataclasses import asdict
from typing import Dict

from studio.app.common.core.utils.config_handler import ConfigWriter
from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.workflow.workflow import Edge, Node, WorkflowConfig
from studio.app.common.core.workflow.workflow_builder import WorkflowConfigBuilder
from studio.app.common.core.workflow.workflow_reader import WorkflowConfigReader
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

    def write(self) -> None:
        config_filepath = join_filepath(
            [
                DIRPATH.OUTPUT_DIR,
                self.workspace_id,
                self.unique_id,
                DIRPATH.WORKFLOW_YML,
            ]
        )

        if os.path.exists(config_filepath):
            config = WorkflowConfigReader.read(config_filepath)
            self.builder.set_config(config)
        else:
            self.create_config()

        ConfigWriter.write(
            dirname=join_filepath(
                [DIRPATH.OUTPUT_DIR, self.workspace_id, self.unique_id]
            ),
            filename=DIRPATH.WORKFLOW_YML,
            config=asdict(self.builder.build()),
        )

    def create_config(self) -> WorkflowConfig:
        return (
            self.builder.set_node_dict(self.nodeDict)
            .set_edge_dict(self.edgeDict)
            .build()
        )
