import pytest
import os
import shutil

from optinist.api.experiment.experiment import ExptConfig, ExptFunction
from optinist.api.experiment.experiment_writer import ExptConfigWriter
from optinist.api.workflow.workflow import Edge, Node, NodeData, RunItem


node_data = NodeData(
    label="a",
    param={},
    path="",
    type=""
)

nodeDict = {
    "test1": Node(
        id="node_id",
        type="a",
        data=node_data,
        position={
            "x": 0,
            "y": 0
        },
        style={
            "border": None,
            "borderRadius": 0,
            "height": 100,
            "padding": 0,
            "width": 180,
        },
    )
}

edgeDict = {
    "test2": Edge(
        id="edge_id",
        type="a",
        animated=False,
        source="",
        sourceHandle="",
        target="",
        targetHandle="",
        style={},
    )
}

dirpath = "/tmp/optinist/output/unique_id"


def test_create_config() -> ExptConfig:
    runItem = RunItem(
        name="New Flow",
        nodeDict=nodeDict,
        edgeDict=edgeDict,
        snakemakeParam={},
        nwbParam={},
        forceRunList=[],
    )

    expt_config = ExptConfigWriter.create_config(
        "test_id",
        runItem.name,
        runItem.nodeDict,
        runItem.edgeDict,
    )

    assert isinstance(expt_config, ExptConfig)
    assert isinstance(expt_config.function, dict)
    assert len(expt_config.function) == 0

    return expt_config


def test_add_run_info():
    expt_config = ExptConfigWriter.add_run_info(
        expt_config=test_create_config(),
        nodeDict=nodeDict,
        edgeDict=None,
    )

    assert len(expt_config.nodeDict) == 1


def test_function_from_nodeDict():
    expt_function = ExptConfigWriter.function_from_nodeDict(nodeDict)

    assert isinstance(expt_function, dict)
    assert isinstance(expt_function["node_id"], ExptFunction)


def test_new_write():
    if os.path.exists(dirpath):
        shutil.rmtree(dirpath)

    ExptConfigWriter.write(
        unique_id="unique_id",
        name="name",
        nodeDict=nodeDict,
        edgeDict=edgeDict,
    )

    assert os.path.exists(f"{dirpath}/experiment.yaml")


def test_write_add():
    ExptConfigWriter.write(
        unique_id="unique_id",
        name="name",
        nodeDict=nodeDict,
        edgeDict=edgeDict,
    )

    assert os.path.exists(f"{dirpath}/experiment.yaml")

    os.remove(f"{dirpath}/experiment.yaml")
