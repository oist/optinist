import os
import shutil
import time

from studio.app.common.core.experiment.experiment import ExptFunction
from studio.app.common.core.rules.runner import Runner
from studio.app.common.core.workflow.workflow import Message, NodeRunStatus
from studio.app.common.core.workflow.workflow_result import (
    NodeResult,
    WorkflowMonitor,
    WorkflowResult,
)
from studio.app.dir_path import DIRPATH

workspace_id = "default"
unique_id = "result_test"
node_id_list = ["func1", "func2"]
node_1st = node_id_list[0]

workflow_dirpath = f"{DIRPATH.DATA_DIR}/output_test/{workspace_id}/{unique_id}"
output_dirpath = f"{DIRPATH.OUTPUT_DIR}/{workspace_id}/{unique_id}"
pickle_path = (
    f"{DIRPATH.DATA_DIR}/output_test/{workspace_id}/{unique_id}/{node_1st}/func1.pkl"
)


def test_WorkflowResult_get_success():
    shutil.copytree(
        workflow_dirpath,
        output_dirpath,
        dirs_exist_ok=True,
    )

    # Write pid_file
    Runner.write_pid_file(
        output_dirpath, "xxxx_dummy_func", "xxxx_dummy_func_script.py"
    )

    output = WorkflowResult(workspace_id=workspace_id, unique_id=unique_id).observe(
        node_id_list
    )

    assert isinstance(output, dict)
    assert output[node_1st].status == "success"


def test_NodeResult_get():
    assert os.path.exists(pickle_path)

    expt_function = ExptFunction(
        unique_id=unique_id,
        name=node_1st,
        success=NodeRunStatus.RUNNING.value,
        hasNWB=False,
    )

    output = NodeResult(
        workspace_id=workspace_id,
        unique_id=unique_id,
        node_id=node_1st,
    ).observe(expt_function)

    assert isinstance(output, Message)


def test_WorkflowResult_get_error():
    shutil.copytree(
        workflow_dirpath,
        output_dirpath,
        dirs_exist_ok=True,
    )

    # Write pid_file file (causes timeout error)
    pid_file_create_time = time.time() - WorkflowMonitor.PROCESS_SNAKEMAKE_WAIT_TIMEOUT
    Runner.write_pid_file(
        output_dirpath,
        "xxxx_dummy_func",
        "xxxx_dummy_func_script.py",
        pid_file_create_time,
    )

    output = WorkflowResult(workspace_id=workspace_id, unique_id=unique_id).observe(
        node_id_list
    )

    assert isinstance(output, dict)
    assert output[node_1st].status == "error"
