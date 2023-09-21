import os
import shutil

from studio.app.common.core.snakemake.smk import ForceRun, SmkParam
from studio.app.common.core.snakemake.snakemake_executor import (
    delete_dependencies,
    snakemake_execute,
)
from studio.app.common.core.workflow.workflow import Edge, Node, NodeData
from studio.app.dir_path import DIRPATH

tiff_filename = "test.tif"

workspace_id = "default"
unique_id = "snakemake"

smk_param = SmkParam(
    use_conda=False,
    cores=2,
    forceall=True,
    forcetargets=True,
    lock=False,
)

shutil.copytree(
    f"{DIRPATH.DATA_DIR}/output_test/{workspace_id}/{unique_id}",
    f"{DIRPATH.DATA_DIR}/output/{workspace_id}/{unique_id}",
    dirs_exist_ok=True,
)


def test_snakemake_execute():
    snakemake_execute(workspace_id, unique_id, smk_param)


nodeDict = {
    "suite2p_file_convert": Node(
        id="suite2p_file_convert",
        type="a",
        data=NodeData(label="suite2p_file_convert", param={}, path="", type=""),
        position={"x": 0, "y": 0},
        style={
            "border": None,
            "borderRadius": 0,
            "height": 100,
            "padding": 0,
            "width": 180,
        },
    ),
    "suite2p_roi": Node(
        id="suite2p_roi",
        type="a",
        data=NodeData(label="suite2p_roi", param={}, path="", type=""),
        position={"x": 0, "y": 0},
        style={
            "border": None,
            "borderRadius": 0,
            "height": 100,
            "padding": 0,
            "width": 180,
        },
    ),
}


edgeDict = {
    "edge1": Edge(
        id="edge1",
        source="input_0",
        animated=False,
        sourceHandle="",
        style={},
        target="suite2p_file_convert",
        targetHandle="",
        type={},
    ),
    "edge2": Edge(
        id="edge2",
        source="suite2p_file_convert",
        animated=False,
        sourceHandle="",
        style={},
        target="suite2p_roi",
        targetHandle="",
        type={},
    ),
}

output_dirpath = f"{DIRPATH.DATA_DIR}/output/default/snakemake"


def test_snakemake_delete_dependencies():
    smk_param.forcerun = [
        ForceRun(
            nodeId="suite2p_roi",
            name="suite2p_roi",
        )
    ]
    delete_dependencies(
        workspace_id=workspace_id,
        unique_id=unique_id,
        smk_params=smk_param,
        nodeDict=nodeDict,
        edgeDict=edgeDict,
    )

    assert not os.path.exists(f"{output_dirpath}/suite2p_roi/suite2p_roi.pkl")


def test_snakemake_delete_dependencies_file_convert():
    test_snakemake_execute()

    smk_param.forcerun = [
        ForceRun(
            nodeId="suite2p_file_convert",
            name="suite2p_file_convert",
        )
    ]
    delete_dependencies(
        workspace_id=workspace_id,
        unique_id=unique_id,
        smk_params=smk_param,
        nodeDict=nodeDict,
        edgeDict=edgeDict,
    )

    assert not os.path.exists(
        f"{output_dirpath}/suite2p_file_convert/suite2p_file_convert.pkl"
    )
    assert not os.path.exists(f"{output_dirpath}/suite2p_roi/suite2p_roi.pkl")


def test_error_snakemake_execute():
    smk_param.use_conda = True
    snakemake_execute(workspace_id, unique_id, smk_param)

    error_log_filepath = f"{DIRPATH.OUTPUT_DIR}/default/snakemake/error.log"

    assert os.path.exists(error_log_filepath)
