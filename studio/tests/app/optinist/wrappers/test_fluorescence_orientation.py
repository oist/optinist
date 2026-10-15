"""FluoData stays (roi, time) in memory. NWBCreater.fluorescence writes the
RoiResponseSeries as (time, roi) by transposing the axis that matches the rois
region, so writers pass their array as is and configs replayed from pickles that
stored (time, roi) still land in the NWB orientation.
"""
import shutil
import sys
import types
import warnings
from pathlib import Path

import numpy as np
import pytest
from pynwb import NWBHDF5IO

from studio.app.common.core.utils.config_handler import ConfigReader
from studio.app.common.core.workflow.workflow_params import read_default_params
from studio.app.common.dataclass import ImageData
from studio.app.dir_path import CORE_PARAM_PATH, DIRPATH
from studio.app.optinist.core.nwb.nwb import NWBDATASET
from studio.app.optinist.core.nwb.nwb_creater import save_nwb
from studio.app.optinist.dataclass import FluoData, IscellData, Suite2pData
from studio.app.optinist.wrappers.lccd import lccd_detection
from studio.app.optinist.wrappers.lccd.lccd_detection import lccd_detect
from studio.app.optinist.wrappers.optinist.neural_population_analysis.correlation import (  # noqa: E501
    correlation,
)
from studio.app.optinist.wrappers.optinist.visualize_utils.vacant_roi import vacant_roi
from studio.app.optinist.wrappers.suite2p.spike_deconv import suite2p_spike_deconv

N_ROI, N_TIME, H, W = 3, 40, 8, 8


@pytest.fixture
def output_dir(tmp_path):
    d = Path(DIRPATH.OUTPUT_DIR) / "default" / tmp_path.name / "func"
    d.mkdir(parents=True)
    yield str(d)
    shutil.rmtree(d.parent)


@pytest.fixture
def images(tmp_path):
    frames = np.random.default_rng(0).random((N_TIME, H, W)).astype(np.float32)
    return ImageData(frames, output_dir=str(tmp_path), file_name="img")


@pytest.fixture
def lccd_stub(monkeypatch):
    class FakeLCCD:
        def __init__(self, params):
            pass

        def apply(self, D):
            roi = np.zeros((D.shape[0] * D.shape[1], N_ROI))
            for i in range(N_ROI):
                roi[i * 4 : (i + 1) * 4, i] = 1
            return roi

    fake = types.ModuleType("lccd")
    fake.LCCD = FakeLCCD
    monkeypatch.setitem(
        sys.modules, f"{lccd_detection.__package__}.lccd_python.lccd", fake
    )


def _series(info):
    (series,) = info["nwbfile"][NWBDATASET.FLUORESCENCE].values()
    return series


def _assert_writer_passes_fluo_as_is(info, fluo_key, n_roi):
    for s in _series(info).values():
        assert s["data"].shape == (n_roi, N_TIME)
        assert len(s["region"]) == n_roi
        if n_roi:  # TimeSeriesData replaces an empty array with zeros((0, 0))
            assert s["data"] is info[fluo_key].data


def _round_trip(tmp_path, images, nwb_config):
    input_config = ConfigReader.read(CORE_PARAM_PATH.nwb.value)
    input_config["image_series"]["external_file"] = images
    path = str(tmp_path / "func.nwb")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        save_nwb(path, input_config, nwb_config)
    assert not [w for w in caught if "oriented incorrectly" in str(w.message)]
    (function_id,) = nwb_config[NWBDATASET.FLUORESCENCE]
    with NWBHDF5IO(path, "r") as io:
        fluo = io.read().processing["ophys"][function_id]
        return {k: s.data[:] for k, s in fluo.roi_response_series.items()}


def test_vacant_roi_passes_fluo_as_is(images, output_dir):
    info = vacant_roi(images, output_dir, params={"f0_frames": 5, "f0_percentile": 8})
    _assert_writer_passes_fluo_as_is(info, "fluorescence", n_roi=0)


def test_lccd_passes_fluo_as_is(images, output_dir, lccd_stub):
    info = lccd_detect(images, output_dir, params={"f0_frames": 5, "f0_percentile": 8})
    _assert_writer_passes_fluo_as_is(info, "fluorescence", n_roi=N_ROI)


def test_spike_deconv_passes_fluo_as_is(output_dir, monkeypatch):
    F = np.random.default_rng(1).random((N_ROI, N_TIME))
    fake = types.ModuleType("suite2p")
    fake.default_ops = lambda: {
        "neucoeff": 0.0,
        "baseline": "maximin",
        "win_baseline": 60.0,
        "sig_baseline": 10.0,
        "prctile_baseline": 8.0,
        "batch_size": 500,
        "tau": 1.0,
    }
    fake.extraction = types.SimpleNamespace(
        preprocess=lambda F, **kw: F, oasis=lambda F, **kw: F
    )
    monkeypatch.setitem(sys.modules, "suite2p", fake)

    ops = Suite2pData({"F": F, "Fneu": np.zeros_like(F), "fs": 30.0, "stat": []})
    info = suite2p_spike_deconv(ops, output_dir, params={})
    _assert_writer_passes_fluo_as_is(info, "spks", n_roi=N_ROI)


def test_nwb_file_is_time_first(images, output_dir, lccd_stub, tmp_path):
    info = lccd_detect(images, output_dir, params={"f0_frames": 5, "f0_percentile": 8})
    fluo = info["fluorescence"].data

    (written,) = _round_trip(tmp_path, images, info["nwbfile"]).values()
    assert written.shape == (N_TIME, N_ROI)
    assert np.array_equal(written, fluo.T)


def test_nwb_file_is_time_first_for_legacy_time_first_config(
    images, output_dir, lccd_stub, tmp_path
):
    info = lccd_detect(images, output_dir, params={"f0_frames": 5, "f0_percentile": 8})
    fluo = info["fluorescence"].data
    _series(info)["Fluorescence"]["data"] = fluo.T.copy()

    (written,) = _round_trip(tmp_path, images, info["nwbfile"]).values()
    assert written.shape == (N_TIME, N_ROI)
    assert np.array_equal(written, fluo.T)


def test_nwb_file_is_time_first_with_zero_rois(images, output_dir, tmp_path):
    info = vacant_roi(images, output_dir, params={"f0_frames": 5, "f0_percentile": 8})

    (written,) = _round_trip(tmp_path, images, info["nwbfile"]).values()
    assert written.shape == (N_TIME, 0)


def test_correlation_default_correlates_rois_not_timepoints(output_dir):
    params = read_default_params("correlation")
    assert params == {"transpose": False}

    fluo = FluoData(np.random.default_rng(2).random((N_ROI, N_TIME)), file_name="f")
    info = correlation(fluo, output_dir, params=params)

    assert info["corr"].data.shape == (N_ROI, N_ROI)


def test_correlation_filters_rois_with_iscell(output_dir):
    data = np.random.default_rng(3).random((5, N_TIME))
    iscell = IscellData(np.array([1, 0, 1, 0, 0]))
    info = correlation(
        FluoData(data, file_name="f"),
        output_dir,
        iscell=iscell,
        params={"transpose": False},
    )

    corr = info["corr"].data
    assert corr.shape == (2, 2)
    assert np.isclose(corr[0, 1], np.corrcoef(data[[0, 2]])[0, 1])
    assert np.isnan(corr[0, 0]) and np.isnan(corr[1, 1])


def test_correlation_handles_a_single_roi(output_dir):
    fluo = FluoData(np.random.default_rng(4).random((1, N_TIME)), file_name="f")
    corr = correlation(fluo, output_dir, params={"transpose": False})["corr"].data
    assert corr.shape == (1, 1)
    assert np.isnan(corr[0, 0])
