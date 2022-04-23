import os
from datetime import datetime
from dataclasses import asdict
from typing import Dict, List

from optinist.api.dir_path import DIRPATH
from optinist.api.experiment.experiment import ExptConfig, ExptFunction
from optinist.api.experiment.experiment_reader import ExptConfigReader
from optinist.api.utils.filepath_creater import join_filepath
from optinist.api.config.config_writer import ConfigWriter
from optinist.api.workflow.workflow import Edge, Node


class ExptConfigWriter:
    @classmethod
    def write(cls, unique_id, name, nodeList, edgeList):
        exp_filepath = join_filepath([DIRPATH.OUTPUT_DIR, unique_id, DIRPATH.EXPERIMENT_YML])
        if os.path.exists(exp_filepath):
            exp_config = ExptConfigReader.read(exp_filepath)
            exp_config = _add_run_info(exp_config, nodeList, edgeList)
        else:
            exp_config = _create(unique_id, name, nodeList, edgeList)

        exp_config.function = _create_function_from_nodeList(nodeList)

        ConfigWriter.write(
            dirname=join_filepath([DIRPATH.OUTPUT_DIR, unique_id]),
            filename=DIRPATH.EXPERIMENT_YML,
            config=asdict(exp_config),
        )


def _create(unique_id, name, nodeList, edgeList) -> ExptConfig:
    return ExptConfig(
        unique_id=unique_id,
        name=name,
        timestamp=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        nodeList=nodeList,
        edgeList=edgeList,
        function={},
    )


def _add_run_info(expt_config: ExptConfig, nodeList: List[Node], edgeList: List[Edge]) -> ExptConfig:
    expt_config.timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # # 関数を追加の可能性
    expt_config.nodeList = nodeList
    expt_config.edgeList = edgeList

    for node in nodeList:
        expt_config.function[node.id] = ExptFunction(
            unique_id=node.id,
            name=node.data.label,
            success="success" if node.data.type == "input" else "running",
        )
    return expt_config


def _create_function_from_nodeList(nodeList: List[Node]) -> Dict[str, ExptFunction]:
    return {
        node.id: ExptFunction(
            unique_id=node.id,
            name=node.data.label,
            success="success" if node.data.type == "input" else "running",
        )
        for node in nodeList
    }
