from studio.app.common.core.workflow.workflow import (
    Edge,
    Node,
    NodeData,
    NodePosition,
    Style,
    WorkflowConfig,
)
from studio.app.common.core.workflow.workflow_reader import WorkflowConfigReader
from studio.app.dir_path import DIRPATH

workflow_filepath = f"{DIRPATH.DATA_DIR}/output_test/0123/workflow.yaml"


def test_read():
    workflow_config = WorkflowConfigReader.read(workflow_filepath)

    assert isinstance(workflow_config, WorkflowConfig)
    assert isinstance(workflow_config.nodeDict, dict)
    assert isinstance(workflow_config.nodeDict["input_0"].id, str)
    assert isinstance(workflow_config.nodeDict["input_0"].type, str)
    assert isinstance(workflow_config.nodeDict["input_0"].data, NodeData)
    assert isinstance(workflow_config.nodeDict["input_0"].position, NodePosition)
    assert isinstance(workflow_config.nodeDict["input_0"].style, Style)

    assert isinstance(workflow_config.edgeDict, dict)

    edgeDictId = (
        "reactflow__edge-input_0input_0--image--ImageData-suite2p_file_"
        "convert_pi2bgrsd6msuite2p_file_convert_pi2bgrsd6m--image--ImageData"
    )
    assert isinstance(workflow_config.edgeDict[edgeDictId].id, str)
    assert isinstance(workflow_config.edgeDict[edgeDictId].type, str)
    assert isinstance(workflow_config.edgeDict[edgeDictId].animated, bool)
    assert isinstance(workflow_config.edgeDict[edgeDictId].source, str)
    assert isinstance(workflow_config.edgeDict[edgeDictId].sourceHandle, str)
    assert isinstance(workflow_config.edgeDict[edgeDictId].target, str)
    assert isinstance(workflow_config.edgeDict[edgeDictId].targetHandle, str)
    assert isinstance(workflow_config.edgeDict[edgeDictId].style, Style)


def test_read_nodeDict():
    node_data = {
        "label": "a",
        "param": {},
        "path": "",
        "type": "",
    }

    position = {
        "x": 0,
        "y": 0,
    }

    nodeDict_config = {
        "sample1": {
            "id": "a",
            "type": "a",
            "data": node_data,
            "position": position,
            "style": {},
        }
    }

    nodeDict = WorkflowConfigReader.read_nodeDict(nodeDict_config)

    assert isinstance(nodeDict, dict)
    assert isinstance(nodeDict["sample1"], Node)


def test_read_edgeDict():
    edgeDict_config = {
        "sample1": {
            "id": "a",
            "type": "a",
            "animated": False,
            "source": "a",
            "sourceHandle": "a",
            "target": "a",
            "targetHandle": "a",
            "style": {},
        }
    }

    edgeDict = WorkflowConfigReader.read_edgeDict(edgeDict_config)

    assert isinstance(edgeDict, dict)
    assert isinstance(edgeDict["sample1"], Edge)
