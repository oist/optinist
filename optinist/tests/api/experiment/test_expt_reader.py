from optinist.api.experiment.experiment import ExptConfig, ExptFunction
from optinist.api.experiment.experiment_reader import ExptConfigReader
from optinist.api.workflow.workflow import Edge, Node, NodeData, NodePosition, Style

expt_filepath = "/tmp/optinist/output/0123/experiment.yaml"


def test_read():
    exp_config = ExptConfigReader.read(expt_filepath)

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

    edgeDictId = "reactflow__edge-input_0input_0--image--ImageData-suite2p_file_ \
        convert_pi2bgrsd6msuite2p_file_convert_pi2bgrsd6m--image--ImageData"
    assert isinstance(exp_config.edgeDict[edgeDictId].id, str)
    assert isinstance(exp_config.edgeDict[edgeDictId].type, str)
    assert isinstance(exp_config.edgeDict[edgeDictId].animated, bool)
    assert isinstance(exp_config.edgeDict[edgeDictId].source, str)
    assert isinstance(exp_config.edgeDict[edgeDictId].sourceHandle, str)
    assert isinstance(exp_config.edgeDict[edgeDictId].target, str)
    assert isinstance(exp_config.edgeDict[edgeDictId].targetHandle, str)
    assert isinstance(exp_config.edgeDict[edgeDictId].style, Style)


def test_read_function():
    func_config = {
        "sample1": {
            "unique_id": "a",
            "name": "a",
            "success": "success",
            "hasNWB": False,
        }
    }

    function = ExptConfigReader.read_function(func_config)

    assert isinstance(function, dict)
    assert isinstance(function["sample1"], ExptFunction)


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

    nodeDict = ExptConfigReader.read_nodeDict(nodeDict_config)

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

    edgeDict = ExptConfigReader.read_edgeDict(edgeDict_config)

    assert isinstance(edgeDict, dict)
    assert isinstance(edgeDict["sample1"], Edge)
