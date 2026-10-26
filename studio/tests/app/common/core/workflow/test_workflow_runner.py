import os
import shutil

import pytest

from studio.app.common.core.experiment.experiment_reader import ExptConfigReader
from studio.app.common.core.workflow.workflow import (
    Edge,
    Node,
    NodeData,
    RunItem,
    WorkflowRunStatus,
)
from studio.app.common.core.workflow.workflow_runner import WorkflowRunner
from studio.app.dir_path import DIRPATH

workspace_id = "default"
unique_id = "workflow_test"

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


dirpath = f"{DIRPATH.OUTPUT_DIR}/{workspace_id}/{unique_id}"


def test_finish_workflow_without_run():
    if os.path.exists(dirpath):
        shutil.rmtree(dirpath)

    runItem = RunItem(
        name="New Flow",
        nodeDict=nodeDict,
        edgeDict=edgeDict,
        snakemakeParam={},
        nwbParam={},
        forceRunList=[],
    )

    WorkflowRunner(workspace_id, unique_id, runItem).finish_workflow_without_run()

    assert os.path.exists(f"{dirpath}/experiment.yaml")
    assert os.path.exists(f"{dirpath}/workflow.yaml")

    exp_config = ExptConfigReader.read(workspace_id, unique_id)

    assert exp_config.success == WorkflowRunStatus.SUCCESS.value


def _stale_run_item():
    # dpca has no top-level "transpose" parameter: it lives under I/O
    stale = Node(
        id="algo",
        type="AlgorithmNode",
        data=NodeData(
            label="dpca",
            param={"transpose": {"type": "child", "value": False}},
            path="",
            type="",
        ),
        position={"x": 0, "y": 0},
        style={},
    )
    return RunItem(
        name="New Flow",
        nodeDict={"algo": stale},
        edgeDict={},
        snakemakeParam={},
        nwbParam={},
        forceRunList=[],
    )


def test_a_stale_parameter_tree_is_refused_before_any_record_is_written():
    shutil.rmtree(dirpath, ignore_errors=True)

    with pytest.raises(KeyError, match="Workflow yaml error"):
        WorkflowRunner(workspace_id, unique_id, _stale_run_item())

    assert not os.path.exists(dirpath)


def test_a_refused_rerun_leaves_the_existing_record_as_it_was():
    shutil.rmtree(dirpath, ignore_errors=True)
    runItem = RunItem(
        name="New Flow",
        nodeDict=nodeDict,
        edgeDict=edgeDict,
        snakemakeParam={},
        nwbParam={},
        forceRunList=[],
    )
    WorkflowRunner(workspace_id, unique_id, runItem).finish_workflow_without_run()
    names = ("experiment.yaml", "workflow.yaml")
    before = {n: open(f"{dirpath}/{n}", "rb").read() for n in names}

    with pytest.raises(KeyError, match="Workflow yaml error"):
        WorkflowRunner(workspace_id, unique_id, _stale_run_item())

    assert {n: open(f"{dirpath}/{n}", "rb").read() for n in names} == before
    assert ExptConfigReader.read(workspace_id, unique_id).success == "success"
    shutil.rmtree(dirpath, ignore_errors=True)
