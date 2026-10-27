import shutil
from pathlib import Path

import numpy as np
import pytest

from studio.app.common.dataclass import TimeSeriesData
from studio.app.dir_path import DIRPATH
from studio.app.optinist.core.nwb.nwb import NWBDATASET
from studio.app.optinist.dataclass import FluoData
from studio.app.optinist.wrappers.optinist.neural_population_analysis.cross_correlation import (  # noqa: E501
    cross_correlation,
)

N_CELL, N_TIME, LAGS = 3, 30, 3
PARAMS = {
    "transpose": False,
    "lags": LAGS,
    "method": "auto",
    "shuffle_sample_number": 2,
    "shuffle_confidence_interval": 0.95,
}


@pytest.fixture
def output_dir(tmp_path):
    d = Path(DIRPATH.OUTPUT_DIR) / "default" / tmp_path.name / "func"
    d.mkdir(parents=True)
    yield str(d)
    shutil.rmtree(d.parent)


def test_every_cell_pair_is_written(output_dir):
    fluo = FluoData(np.random.default_rng(0).random((N_CELL, N_TIME)), file_name="f")
    info = cross_correlation(fluo, output_dir, params=PARAMS)

    pairs = ["0-1", "0-2", "1-2"]
    assert sorted(k for k in info if k != "nwbfile") == sorted(
        pairs + [f"shuffle {p}" for p in pairs]
    )
    n_lags = 2 * LAGS + 1
    for p in pairs:
        assert isinstance(info[p], TimeSeriesData)
        assert info[p].data.shape == (2, n_lags)
        assert info[f"shuffle {p}"].data.shape == (4, n_lags)
        assert info[p].data[0].tolist() == list(range(-LAGS, LAGS + 1))
    (post,) = info["nwbfile"][NWBDATASET.POSTPROCESS].values()
    assert post["mat"].shape == (N_CELL, N_CELL, n_lags)


def test_cross_correlation_leaves_its_input_unchanged(output_dir):
    data = np.random.default_rng(1).random((N_CELL, N_TIME))
    fluo = FluoData(data, file_name="f")
    nwb_fluorescence = {"Fluorescence": {"data": fluo.data}}
    original = data.copy()

    info = cross_correlation(fluo, output_dir, params=PARAMS)

    assert np.array_equal(fluo.data, original)
    assert np.array_equal(nwb_fluorescence["Fluorescence"]["data"], original)
    (post,) = info["nwbfile"][NWBDATASET.POSTPROCESS].values()
    lag0 = original @ original.T
    assert np.allclose(post["mat"][:, :, LAGS], lag0)


def test_cross_correlation_needs_two_cells(output_dir):
    fluo = FluoData(np.random.default_rng(5).random((1, N_TIME)), file_name="f")
    with pytest.raises(ValueError, match="at least 2 cells"):
        cross_correlation(fluo, output_dir, params={"transpose": False})
