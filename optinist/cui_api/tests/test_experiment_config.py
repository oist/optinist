import pytest

from typing import List
from optinist.cui_api.dir_path import DIRPATH
from optinist.cui_api.filepath_creater import (
    create_filepath,
    join_filepath
)
from optinist.cui_api.experiment_config import (
    ExpConfig,
    ExpConfigReader,
    NodeData,
    NodePosition,
    Style,
    create_exp_config
)


from dataclasses import dataclass

@dataclass
class RunItem:
    name: str
    nodeList: list
    edgeList: list
    snakemakeParam: dict
    nwbParam: dict

def test_filepath():
    exp_filepath = create_filepath(
        join_filepath([DIRPATH.ROOT_DIR, "test_data"]),
        DIRPATH.EXPERIMENT_YML
    )
    return exp_filepath


def test_exp_config_reader():
    exp_filepath = test_filepath()
    exp_config = ExpConfigReader.read(exp_filepath)

    assert isinstance(exp_config, ExpConfig)
    assert isinstance(exp_config.timestamp, str)
    assert isinstance(exp_config.name, str)
    assert isinstance(exp_config.unique_id, str)

    assert isinstance(exp_config.function, dict)
    assert isinstance(exp_config.function["0"].unique_id, str)
    assert isinstance(exp_config.function["0"].name, str)
    assert isinstance(exp_config.function["0"].success, str)

    assert isinstance(exp_config.nodeList, list)
    assert isinstance(exp_config.nodeList[0].id, str)
    assert isinstance(exp_config.nodeList[0].type, str)
    assert isinstance(exp_config.nodeList[0].data, NodeData)
    assert isinstance(exp_config.nodeList[0].position, NodePosition)
    assert isinstance(exp_config.nodeList[0].style, Style)

    assert isinstance(exp_config.edgeList, list)
    assert isinstance(exp_config.edgeList[0].id, str)
    assert isinstance(exp_config.edgeList[0].type, str)
    assert isinstance(exp_config.edgeList[0].animated, bool)
    assert isinstance(exp_config.edgeList[0].source, str)
    assert isinstance(exp_config.edgeList[0].sourceHandle, str)
    assert isinstance(exp_config.edgeList[0].target, str)
    assert isinstance(exp_config.edgeList[0].targetHandle, str)
    assert isinstance(exp_config.edgeList[0].style, Style)
    

def test_create_exp_config():
    exp_filepath = test_filepath()

    runItem = RunItem(
        name="New Flow",
        nodeList=[],
        edgeList=[],
        snakemakeParam={},
        nwbParam={},
    )

    exp_config = create_exp_config(
        exp_filepath,
        runItem.name,
        runItem.nodeList,
        runItem.edgeList,
    )

    assert isinstance(exp_config, ExpConfig)
    assert isinstance(exp_config.function, dict)
    assert len(exp_config.function) == 0
