import os
import shutil

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
