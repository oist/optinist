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
    assert isinstance(exp_config.function["input_0"].unique_id, str)
    assert isinstance(exp_config.function["input_0"].name, str)
    assert isinstance(exp_config.function["input_0"].success, str)

    assert isinstance(exp_config.nodeDict, dict)
    assert isinstance(exp_config.nodeDict["input_0"].id, str)
    assert isinstance(exp_config.nodeDict["input_0"].type, str)
    assert isinstance(exp_config.nodeDict["input_0"].data, NodeData)
    assert isinstance(exp_config.nodeDict["input_0"].position, NodePosition)
    assert isinstance(exp_config.nodeDict["input_0"].style, Style)

    assert isinstance(exp_config.edgeDict, dict)
    
    edgeDictId = "reactflow__edge-input_0input_0--image--ImageData-suite2p_file_convert_pi2bgrsd6msuite2p_file_convert_pi2bgrsd6m--image--ImageData"
    assert isinstance(exp_config.edgeDict[edgeDictId].id, str)
    assert isinstance(exp_config.edgeDict[edgeDictId].type, str)
    assert isinstance(exp_config.edgeDict[edgeDictId].animated, bool)
    assert isinstance(exp_config.edgeDict[edgeDictId].source, str)
    assert isinstance(exp_config.edgeDict[edgeDictId].sourceHandle, str)
    assert isinstance(exp_config.edgeDict[edgeDictId].target, str)
    assert isinstance(exp_config.edgeDict[edgeDictId].targetHandle, str)
    assert isinstance(exp_config.edgeDict[edgeDictId].style, Style)
