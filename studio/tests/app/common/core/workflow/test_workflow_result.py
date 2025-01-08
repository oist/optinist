import os
import shutil

from studio.app.common.core.rules.runner import Runner
from studio.app.common.core.workflow.workflow import Message
from studio.app.common.core.workflow.workflow_result import NodeResult, WorkflowResult
from studio.app.dir_path import DIRPATH

workspace_id = "default"
unique_id = "result_test"
node_id_list = ["func1", "func2"]

workflow_dirpath = f"{DIRPATH.DATA_DIR}/output_test/{workspace_id}/{unique_id}"
output_dirpath = f"{DIRPATH.OUTPUT_DIR}/{workspace_id}/{unique_id}"
pickle_path = (
    f"{DIRPATH.DATA_DIR}/output_test/{workspace_id}/{unique_id}/func1/func1.pkl"
)


def test_WorkflowResult_get():
    shutil.copytree(
        workflow_dirpath,
        output_dirpath,
        dirs_exist_ok=True,
    )

    # first, write pid_file
    Runner.write_pid_file(output_dirpath, "xxxx_dummy_func_script.py")

    output = WorkflowResult(workspace_id=workspace_id, unique_id=unique_id).observe(
        node_id_list
    )

    assert isinstance(output, dict)
    assert len(output) == 1


def test_NodeResult_get():
    assert os.path.exists(pickle_path)
    output = NodeResult(
        workspace_id=workspace_id,
        unique_id=unique_id,
        node_id="func1",
        pickle_filepath=pickle_path,
    ).observe()

    assert isinstance(output, Message)
