import os
from datetime import datetime
from dataclasses import asdict

from optinist.api.dir_path import DIRPATH
from optinist.api.experiment.experiment import ExpConfig, ExpFunction
from optinist.api.experiment.experiment_config_reader import ExpConfigReader
from optinist.api.utils.filepath_creater import join_filepath
from optinist.api.config.config_writer import ConfigWriter

from typing import List

from optinist.api.workflow.workflow import Edge, Node


def create_exp_config(unique_id, name, nodeList, edgeList):
    return ExpConfig(
        unique_id=unique_id,
        name=name,
        timestamp=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        nodeList=nodeList,
        edgeList=edgeList,
        function={},
    )


def add_run_info(exp_config: ExpConfig, nodeList: List[Node], edgeList: List[Edge]):
    exp_config.timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # # 関数を追加の可能性
    exp_config.nodeList = nodeList
    exp_config.edgeList = edgeList

    for node in nodeList:
        exp_config.function[node.id] = ExpFunction(
            unique_id=node.id,
            name=node.data.label,
            success="success" if node.data.type == "input" else "running",
        )
    return exp_config


def create_function_from_nodeList(nodeList: List[Node]):
    return {
        node.id: ExpFunction(
            unique_id=node.id,
            name=node.data.label,
            success="success" if node.data.type == "input" else "running",
        )
        for node in nodeList
    }


def exp_config_writer(unique_id, name, nodeList, edgeList):
    exp_filepath = join_filepath([DIRPATH.BASE_DIR, unique_id, DIRPATH.EXPERIMENT_YML])
    if os.path.exists(exp_filepath):
        exp_config = ExpConfigReader.read(exp_filepath)
        exp_config = add_run_info(exp_config, nodeList, edgeList)
    else:
        exp_config = create_exp_config(unique_id, name, nodeList, edgeList)

    exp_config.function = create_function_from_nodeList(nodeList)

    ConfigWriter.write(
        dirname=join_filepath([DIRPATH.BASE_DIR, unique_id]),
        filename=DIRPATH.EXPERIMENT_YML,
        config=asdict(exp_config),
    )
