import shutil

import pytest
from pynwb import NWBHDF5IO

from studio.app.common.core.mode import MODE
from studio.app.common.core.snakemake.smk import SmkParam
from studio.app.common.core.snakemake.snakemake_executor import snakemake_execute
from studio.app.common.core.utils.pickle_handler import PickleReader
from studio.app.dir_path import DIRPATH

workspace_id = "default"
unique_id = "smk_exec_lccd_dpca"
output_dirpath = f"{DIRPATH.OUTPUT_DIR}/{workspace_id}/{unique_id}"
HEATMAP_NAMES = [f"{f}-component{c}" for f in ["t", "b", "c", "tbc"] for c in [0, 1]]


@pytest.mark.lighter_processing
def test_lccd_dpca_workflow_on_real_data(client):
    """lccd on the test movie, then dPCA with the default params and real behavior."""
    shutil.copytree(
        f"{DIRPATH.DATA_DIR}/output_test/{workspace_id}/{unique_id}",
        output_dirpath,
        dirs_exist_ok=True,
    )
    MODE.reset_mode(is_standalone=True)

    assert snakemake_execute(
        workspace_id, unique_id, SmkParam(use_conda=True, cores=2, forceall=True)
    )

    result = PickleReader.read(f"{output_dirpath}/dpca_1234/dpca.pkl")
    assert PickleReader.check_is_valid_node_pickle(result), result
    assert [k for k in result if k != "nwbfile"] == HEATMAP_NAMES
    heat = result["t-component0"]
    assert heat.data.shape == (4, 20)
    assert list(heat.columns) == list(range(-10, 10))

    with NWBHDF5IO(f"{output_dirpath}/whole.nwb", "r") as io:
        optinist = io.read().processing["optinist"]
        assert optinist["dpca_1234_tbc"].data.shape == (8, 20, 2, 2)
        assert optinist["dpca_1234_trial_counts"].data[:].sum() > 0
