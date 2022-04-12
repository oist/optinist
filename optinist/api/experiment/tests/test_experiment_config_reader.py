import pytest

from optinist.api.dir_path import DIRPATH
from optinist.api.experiment.experiment_reader import ExptConfigReader
from optinist.api.utils.filepath_creater import (
    create_filepath,
    join_filepath
)
from optinist.api.experiment.experiment import ExptConfig
from optinist.api.workflow.workflow import NodeData, NodePosition, Style


def test_filepath():
    exp_filepath = create_filepath(
        join_filepath([DIRPATH.ROOT_DIR, "test_data"]),
        DIRPATH.EXPERIMENT_YML
    )
    return exp_filepath


def test_read():
    exp_filepath = test_filepath()
    exp_config = ExptConfigReader.read(exp_filepath)

    assert isinstance(exp_config, ExptConfig)
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
