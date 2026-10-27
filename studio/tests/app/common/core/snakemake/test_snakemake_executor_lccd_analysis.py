import itertools
import shutil

import numpy as np
import pytest
from pynwb import NWBHDF5IO

from studio.app.common.core.mode import MODE
from studio.app.common.core.snakemake.smk import SmkParam
from studio.app.common.core.snakemake.snakemake_executor import snakemake_execute
from studio.app.common.core.utils.pickle_handler import PickleReader
from studio.app.dir_path import DIRPATH

workspace_id = "default"
unique_id = "smk_exec_lccd_analysis"
output_dirpath = f"{DIRPATH.OUTPUT_DIR}/{workspace_id}/{unique_id}"


def _node_result(node_id, label):
    result = PickleReader.read(f"{output_dirpath}/{node_id}/{label}.pkl")
    assert PickleReader.check_is_valid_node_pickle(result), f"{node_id}: {result}"
    return result


@pytest.mark.lighter_processing
def test_lccd_analysis_workflow_on_real_data(client):
    """lccd on the test movie, then correlation and cross_correlation, end to end."""
    shutil.copytree(
        f"{DIRPATH.DATA_DIR}/output_test/{workspace_id}/{unique_id}",
        output_dirpath,
        dirs_exist_ok=True,
    )
    MODE.reset_mode(is_standalone=True)

    assert snakemake_execute(
        workspace_id, unique_id, SmkParam(use_conda=True, cores=2, forceall=True)
    )

    lccd = _node_result("lccd_cell_detection_1234", "lccd_cell_detection")
    fluo = lccd["fluorescence"].data
    cells = np.where(lccd["iscell"].data > 0)[0]
    n_roi, n_time = fluo.shape
    assert n_time == 1000
    assert len(cells) >= 2, "test movie should yield at least two cells"

    corr = _node_result("correlation_1234", "correlation")["corr"].data
    assert corr.shape == (len(cells), len(cells))
    assert np.isclose(corr[0, 1], np.corrcoef(fluo[cells[:2]])[0, 1])

    cross = _node_result("cross_correlation_1234", "cross_correlation")
    pairs = [f"{i}-{j}" for i, j in itertools.combinations(range(len(cells)), 2)]
    assert sorted(k for k in cross if k != "nwbfile") == sorted(
        pairs + [f"shuffle {p}" for p in pairs]
    )

    with NWBHDF5IO(f"{output_dirpath}/whole.nwb", "r") as io:
        ophys = io.read().processing["ophys"]
        series = ophys["lccd_cell_detection_1234"].roi_response_series["Fluorescence"]
        assert series.data.shape == (n_time, n_roi)
        assert np.allclose(series.data[:], fluo.T)
