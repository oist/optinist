import pytest
import os
from optinist.api.workflow.workflow import Message

from optinist.api.workflow.workflow_result import NodeResult, WorkflowResult


unique_id = "result_test"
node_id_list = ["func1", "func2"]

workflow_dirpath = f"/tmp/optinist/output/{unique_id}"
pickle_path = f"/tmp/optinist/output/{unique_id}/func1/func1.pkl"


def test_WorkflowResult_get():
    output = WorkflowResult(
        unique_id=unique_id
    ).get(node_id_list)

    assert isinstance(output, dict)
    assert len(output) == 1


def test_NodeResult_get():
    assert os.path.exists(pickle_path)
    output = NodeResult(
        workflow_dirpath=workflow_dirpath,
        node_id="func1",
        pickle_filepath=pickle_path,
    ).get()

    assert isinstance(output, Message)
