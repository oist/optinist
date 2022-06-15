import pytest
import shutil

from optinist.api.workflow.workflow_result import WorkflowResult


def test_WorkflowResult_get():
    output = WorkflowResult("test_data").get(["test_workflow"])

    assert isinstance(output, dict)
    assert len(output) == 0

    # shutil.copytree("./test1", "./test2")
    # output = WorkflowResult("test_data").get(["test_workflow"])

    # assert isinstance(output, dict)
    # assert len(output) == 0
