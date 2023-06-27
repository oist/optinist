import os
import shutil

from optinist.services.dir_path import DIRPATH
from optinist.services.experiment.experiment import ExptConfig, ExptFunction
from optinist.services.experiment.experiment_writer import ExptConfigWriter
from optinist.services.workflow.workflow import Edge, Node, NodeData, RunItem

node_data = NodeData(label="a", param={}, path="", type="")

nodeDict = {
    "test1": Node(
        id="node_id",
        type="a",
        data=node_data,
        position={"x": 0, "y": 0},
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

dirpath = f"{DIRPATH.OPTINIST_DIR}/output/unique_id"


def test_create_config() -> ExptConfig:
    runItem = RunItem(
        name="New Flow",
        nodeDict=nodeDict,
        edgeDict=edgeDict,
        snakemakeParam={},
        nwbParam={},
        forceRunList=[],
    )

    expt_config = ExptConfigWriter(
        unique_id="test_id",
        name=runItem.name,
        nodeDict=runItem.nodeDict,
        edgeDict=runItem.edgeDict,
    ).create_config()

    assert isinstance(expt_config, ExptConfig)
    assert isinstance(expt_config.function, dict)
    assert len(expt_config.function) == 0

    assert expt_config


def test_add_run_info():
    expt_config = ExptConfigWriter(
        unique_id="",
        name="",
        nodeDict=nodeDict,
        edgeDict=None,
    ).add_run_info()

    assert len(expt_config.nodeDict) == 1


def test_function_from_nodeDict():
    expt_config = ExptConfigWriter(
        unique_id="",
        name="",
        nodeDict=nodeDict,
        edgeDict=edgeDict,
    ).function_from_nodeDict()

    assert isinstance(expt_config.function, dict)
    assert isinstance(expt_config.function["node_id"], ExptFunction)


def test_new_write():
    if os.path.exists(dirpath):
        shutil.rmtree(dirpath)

    ExptConfigWriter(
        unique_id="unique_id",
        name="name",
        nodeDict=nodeDict,
        edgeDict=edgeDict,
    ).write()

    assert os.path.exists(f"{dirpath}/experiment.yaml")


def test_write_add():
    ExptConfigWriter(
        unique_id="unique_id",
        name="name",
        nodeDict=nodeDict,
        edgeDict=edgeDict,
    ).write()

    assert os.path.exists(f"{dirpath}/experiment.yaml")

    os.remove(f"{dirpath}/experiment.yaml")
