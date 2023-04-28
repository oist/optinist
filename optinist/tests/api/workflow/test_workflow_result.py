import os
import shutil

from optinist.api.dir_path import DIRPATH
from optinist.api.workflow.workflow import Message
from optinist.api.workflow.workflow_result import NodeResult, WorkflowResult

unique_id = "result_test"
node_id_list = ["func1", "func2"]

workflow_dirpath = f"{DIRPATH.OPTINIST_DIR}/output_test/{unique_id}"
pickle_path = f"{DIRPATH.OPTINIST_DIR}/output_test/{unique_id}/func1/func1.pkl"


def test_WorkflowResult_get():
    shutil.copytree(
        workflow_dirpath,
        f"{DIRPATH.OPTINIST_DIR}/output/{unique_id}",
        dirs_exist_ok=True,
    )
    output = WorkflowResult(unique_id=unique_id).get(node_id_list)

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
