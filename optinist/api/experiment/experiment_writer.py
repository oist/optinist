import os
from dataclasses import asdict
from datetime import datetime
from typing import Dict

from optinist.api.config.config_writer import ConfigWriter
from optinist.api.dir_path import DIRPATH
from optinist.api.experiment.experiment import ExptConfig, ExptFunction
from optinist.api.experiment.experiment_builder import ExptConfigBuilder
from optinist.api.experiment.experiment_reader import ExptConfigReader
from optinist.api.utils.filepath_creater import join_filepath
from optinist.api.workflow.workflow import Edge, Node


class ExptConfigWriter:
    def __init__(
        self,
        unique_id: str,
        name: str,
        nodeDict: Dict[str, Node],
        edgeDict: Dict[str, Edge],
    ) -> None:
        self.unique_id = unique_id
        self.name = name
        self.nodeDict = nodeDict
        self.edgeDict = edgeDict

        self.builder = ExptConfigBuilder()

    def write(self) -> None:
        expt_filepath = join_filepath(
            [DIRPATH.OUTPUT_DIR, self.unique_id, DIRPATH.EXPERIMENT_YML]
        )
        if os.path.exists(expt_filepath):
            expt_config = ExptConfigReader.read(expt_filepath)
            self.builder.set_config(expt_config)
            self.add_run_info()
        else:
            self.create_config()

        self.function_from_nodeDict()

        ConfigWriter.write(
            dirname=join_filepath([DIRPATH.OUTPUT_DIR, self.unique_id]),
            filename=DIRPATH.EXPERIMENT_YML,
            config=asdict(self.builder.build()),
        )

    def create_config(self) -> ExptConfig:
        return (
            self.builder.set_unique_id(self.unique_id)
            .set_name(self.name)
            .set_timestamp(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            .set_nodeDict(self.nodeDict)
            .set_edgeDict(self.edgeDict)
            .build()
        )

    def add_run_info(self) -> ExptConfig:
        return (
            self.builder.set_timestamp(
                datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            )  # 時間を更新
            .set_nodeDict(self.nodeDict)
            .set_edgeDict(self.edgeDict)
            .build()
        )

    def function_from_nodeDict(self) -> ExptConfig:
        func_dict = {
            node.id: ExptFunction(
                unique_id=node.id,
                name=node.data.label,
                success="success" if node.data.type == "input" else "running",
                hasNWB=False,
            )
            for node in self.nodeDict.values()
        }
        return self.builder.set_function(func_dict).build()
