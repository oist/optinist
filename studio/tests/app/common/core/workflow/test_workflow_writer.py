import os
import shutil

from studio.app.common.core.workflow.workflow import (
    Edge,
    Node,
    NodeData,
    WorkflowConfig,
)
from studio.app.common.core.workflow.workflow_writer import WorkflowConfigWriter
from studio.app.dir_path import DIRPATH

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


def test_create_config() -> WorkflowConfig:
    workflow_config = WorkflowConfigWriter(
        workspace_id="test_workspace_id",
        unique_id="test_id",
        nodeDict=nodeDict,
        edgeDict=edgeDict,
    ).create_config()

    assert isinstance(workflow_config, WorkflowConfig)


dirpath = f"{DIRPATH.DATA_DIR}/output/workspace_id/unique_id"


def test_write():
    if os.path.exists(dirpath):
        shutil.rmtree(dirpath)

    WorkflowConfigWriter(
        workspace_id="workspace_id",
        unique_id="unique_id",
        nodeDict=nodeDict,
        edgeDict=edgeDict,
    ).write()

    assert os.path.exists(f"{dirpath}/workflow.yaml")
