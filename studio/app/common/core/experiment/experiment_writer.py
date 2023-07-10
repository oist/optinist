import os
from dataclasses import asdict
from datetime import datetime
from typing import Dict

from studio.app.common.core.experiment.experiment import ExptConfig, ExptFunction
from studio.app.common.core.experiment.experiment_builder import ExptConfigBuilder
from studio.app.common.core.experiment.experiment_reader import ExptConfigReader
from studio.app.common.core.utils.config_handler import ConfigWriter
from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.workflow.workflow import Edge, Node
from studio.app.const import DATE_FORMAT
from studio.app.dir_path import DIRPATH


class ExptConfigWriter:
    def __init__(
        self,
        workspace_id: str,
        unique_id: str,
        name: str,
        nodeDict: Dict[str, Node],
        edgeDict: Dict[str, Edge],
    ) -> None:
        self.workspace_id = workspace_id
        self.unique_id = unique_id
        self.name = name
        self.nodeDict = nodeDict
        self.edgeDict = edgeDict

        self.builder = ExptConfigBuilder()

    def write(self) -> None:
        expt_filepath = join_filepath(
            [
                DIRPATH.OUTPUT_DIR,
                self.workspace_id,
                self.unique_id,
                DIRPATH.EXPERIMENT_YML,
            ]
        )
        if os.path.exists(expt_filepath):
            expt_config = ExptConfigReader.read(expt_filepath)
            self.builder.set_config(expt_config)
            self.add_run_info()
        else:
            self.create_config()

        self.function_from_nodeDict()

        ConfigWriter.write(
            dirname=join_filepath(
                [DIRPATH.OUTPUT_DIR, self.workspace_id, self.unique_id]
            ),
            filename=DIRPATH.EXPERIMENT_YML,
            config=asdict(self.builder.build()),
        )

    def create_config(self) -> ExptConfig:
        return (
            self.builder.set_workspace_id(self.workspace_id)
            .set_unique_id(self.unique_id)
            .set_name(self.name)
            .set_started_at(datetime.now().strftime(DATE_FORMAT))
            .set_success("running")
            .set_nodeDict(self.nodeDict)
            .set_edgeDict(self.edgeDict)
            .build()
        )

    def add_run_info(self) -> ExptConfig:
        return (
            self.builder.set_started_at(datetime.now().strftime(DATE_FORMAT))  # 時間を更新
            .set_success("running")
            .set_nodeDict(self.nodeDict)
            .set_edgeDict(self.edgeDict)
            .build()
        )

    def function_from_nodeDict(self) -> ExptConfig:
        func_dict: Dict[str, ExptFunction] = {}
        for node in self.nodeDict.values():
            func_dict[node.id] = ExptFunction(
                unique_id=node.id, name=node.data.label, hasNWB=False, success="running"
            )
            if node.data.type == "input":
                timestamp = datetime.now().strftime(DATE_FORMAT)
                func_dict[node.id].started_at = timestamp
                func_dict[node.id].finished_at = timestamp
                func_dict[node.id].success = "success"

        return self.builder.set_function(func_dict).build()
