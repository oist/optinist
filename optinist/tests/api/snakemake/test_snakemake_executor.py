import os
import shutil

import pytest

from optinist.api.dir_path import DIRPATH
from optinist.api.snakemake.smk import ForceRun, SmkParam
from optinist.api.snakemake.snakemake_executor import (
    delete_dependencies,
    snakemake_execute,
)
from optinist.api.workflow.workflow import Edge, Node, NodeData

tiff_filename = "test.tif"

unique_id = "snakemake"

smk_param = SmkParam(
    use_conda=False,
    cores=2,
    forceall=True,
    forcetargets=True,
    lock=False,
)

shutil.copyfile(
    "/tmp/optinist/config.yaml",
    f"{DIRPATH.ROOT_DIR}/config.yaml",
)


def test_snakemake_execute():
    snakemake_execute(unique_id, smk_param)


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

output_dirpath = "/tmp/optinist/output/snakemake"


def test_snakemake_delete_dependencies():
    smk_param.forcerun = [
        ForceRun(
            nodeId="suite2p_roi",
            name="suite2p_roi",
        )
    ]
    delete_dependencies(
        unique_id=unique_id,
        smk_params=smk_param,
        nodeDict=nodeDict,
        edgeDict=edgeDict,
    )

    assert not os.path.exists(f"{output_dirpath}/suite2p_roi/suite2p_roi.pkl")


def test_snakemake_delete_dependencies():
    test_snakemake_execute()

    smk_param.forcerun = [
        ForceRun(
            nodeId="suite2p_file_convert",
            name="suite2p_file_convert",
        )
    ]
    delete_dependencies(
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
    snakemake_execute(unique_id, smk_param)

    error_log_filepath = "/tmp/optinist/output/snakemake/error.log"

    assert os.path.exists(error_log_filepath)

    os.remove(error_log_filepath)
