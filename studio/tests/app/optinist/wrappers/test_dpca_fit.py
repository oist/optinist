import importlib.util
import inspect
import itertools
import re
import sys
import types
from pathlib import Path

import numpy as np
import pytest
import yaml

from studio.app.dir_path import DIRPATH
from studio.app.optinist.core.nwb.nwb import NWBDATASET
from studio.app.optinist.dataclass import BehaviorData, FluoData

mod = importlib.import_module(
    "studio.app.optinist.wrappers.optinist.dimension_reduction.dpca_fit"
)

REPO_ROOT = Path(__file__).resolve().parents[5]
SAMPLE_CSV = REPO_ROOT / "sample_data/tutorial/input/sample_mouse2p_behavior.csv"
N_ROI, N_TIME = 30, 2000
OUTPUT_DIR = f"{DIRPATH.OUTPUT_DIR}/default/uid/dpca_x"
HEATMAP_NAMES = [f"{f}-component{c}" for f in ["t", "b", "c", "tbc"] for c in [0, 1]]


@pytest.fixture(autouse=True)
def no_sklearn_needed(monkeypatch):
    # The poetry image has no scikit-learn; the guards under test do not need
    # the real standardization, so the local lane runs it as the identity.
    if importlib.util.find_spec("sklearn") is None:
        monkeypatch.setattr(mod, "standard_norm", lambda X, mean, std: X)


@pytest.fixture
def defaults():
    path = Path(inspect.getfile(mod)).parent / "params/dpca.yaml"
    nested = yaml.safe_load(path.read_text())
    flat = {}
    mod.recursive_flatten_params(nested, flat)
    return flat


@pytest.fixture
def behavior():
    return np.loadtxt(SAMPLE_CSV, delimiter=",")


@pytest.fixture
def fluo():
    return np.random.default_rng(0).random((N_ROI, N_TIME))


def _prepare(fluo, behavior, defaults, **over):
    params = {**defaults, **over}
    return mod.prepare_inputs(fluo, behavior, None, params)


def test_build_trials_groups_all_features_and_averages_every_trial():
    rng = np.random.default_rng(1)
    D = rng.random((200, 5))
    triggers = np.arange(20, 180, 10)
    f1 = np.array([0, 1] * 8)
    f2 = np.array([0, 0, 1, 1] * 4)
    f3 = np.array([0] * 12 + [1] * 4)
    mean, counts, levels = mod.build_trials(
        D, triggers, [f1, f2, f3], [0, 1, 2], [-2, 2]
    )

    assert mean.shape == (5, 4, 2, 2, 2)
    assert counts.tolist() == [[[3, 1], [3, 1]], [[3, 1], [3, 1]]]
    assert [list(u) for u in levels] == [[0, 1], [0, 1], [0, 1]]
    cond000 = [t for t, a, b, c in zip(triggers, f1, f2, f3) if a == b == c == 0]
    expected = np.mean([D[t - 2 : t + 2].T for t in cond000], axis=0)
    assert np.allclose(mean[:, :, 0, 0, 0], expected)
    assert np.isfinite(mean).all()


def test_build_trials_drops_windows_outside_recording():
    D = np.random.default_rng(2).random((100, 3))
    triggers = np.array([2, 50, 60, 99])
    f = np.array([0, 1, 0, 1])
    _, counts, _ = mod.build_trials(D, triggers, [f], [0], [-5, 5])
    assert counts.tolist() == [1, 1]
    with pytest.raises(ValueError, match="no trigger window"):
        mod.build_trials(D, np.array([2, 99]), [np.array([0, 1])], [0], [-5, 5])


def test_build_trials_rejects_missing_combination_and_bad_features():
    D = np.zeros((100, 2))
    triggers = np.arange(10, 90, 10)
    with pytest.raises(ValueError, match=r"no trials for .*\(1, 1\)"):
        mod.build_trials(
            D,
            triggers,
            [np.array([0, 0, 0, 0, 1, 1, 1, 1]), np.array([0, 1, 0, 1, 0, 0, 0, 0])],
            [0, 2],
            [-1, 1],
        )
    with pytest.raises(ValueError, match="single value"):
        mod.build_trials(D, triggers, [np.zeros(8)], [4], [-1, 1])
    with pytest.raises(ValueError, match="continuous"):
        mod.build_trials(D, np.arange(2, 98, 4), [np.arange(24) / 23], [3], [-1, 1])


def test_build_trials_caps_the_missing_combination_list():
    D = np.zeros((100, 2))
    triggers = np.arange(10, 70, 10)
    diagonal = np.arange(6)
    with pytest.raises(ValueError) as err:
        mod.build_trials(D, triggers, [diagonal, diagonal], [0, 1], [-1, 1])
    message = str(err.value)
    assert len(re.findall(r"\(\d, \d\)", message)) == mod.MAX_LISTED_COMBINATIONS
    assert "and 20 more" in message


def test_build_trials_warns_on_unbalanced_trial_counts(caplog):
    D = np.zeros((200, 2))
    triggers = np.arange(10, 190, 10)
    f = np.array([0] * 17 + [1])
    with caplog.at_level("WARNING"):
        mod.build_trials(D, triggers, [f], [0], [-1, 1])
    assert "unbalanced trial counts per condition [17, 1]" in caplog.text


def test_calc_trigger_counts_signal_already_high_at_frame_zero():
    signal = np.array([1, 1, 0, 0, 1, 0])
    assert mod.calc_trigger(signal, "up", 0.5).tolist() == [0, 4]
    assert mod.calc_trigger(signal, "cross", 0.5).tolist() == [0, 2, 4, 5]
    assert mod.calc_trigger(signal, "down", 0.5).tolist() == [2, 5]


def test_defaults_run_on_sample_behavior(fluo, behavior, defaults):
    prepared = _prepare(fluo, behavior, defaults)
    assert prepared["mean"].shape == (N_ROI, 20, 2, 2)
    assert prepared["counts"].sum() == 99
    assert prepared["labels"] == "tbc"
    assert prepared["regularizer"] == 0.0
    assert prepared["seed"] == 0
    assert prepared["figure_features"] == ["t", "b", "c", "tbc"]


@pytest.mark.parametrize(
    "over, match",
    [
        ({"transpose": False}, "share the time axis"),
        ({"feature_columns": [3, 9]}, "out of range"),
        ({"trigger_duration": [10, -10]}, "before < after"),
        ({"labels": "tb"}, "labels must be 3"),
        ({"figure_features": ["t", "x"]}, "not marginalizations"),
        ({"regularizer": "auto"}, "regularizer must be a number"),
        ({"regularizer": -0.5}, "regularizer must be a number"),
        ({"n_components": 40}, "exceeds the number of cells"),
        ({"figure_components": [0, 9]}, "figure_components"),
        ({"trigger_threshold": 5}, "no 'up' trigger"),
        ({"feature_columns": [], "labels": "t"}, "at least one behavior column"),
        ({"feature_columns": [0, 7]}, r"no trials for .*\(1.0, -1.0\)"),
        ({"feature_columns": [3, 7]}, "continuous"),
        ({"regularizer": "inf"}, "regularizer must be a number"),
        ({"regularizer": "nan"}, "regularizer must be a number"),
        ({"seed": -1}, "seed must be between"),
        ({"feature_columns": [0, 2, 3], "labels": "tbcd"}, "at most 2 columns"),
        ({"feature_columns": [0.9, 2.1]}, "feature_columns must be whole numbers"),
        ({"trigger_column": 1.9}, "trigger_column must be whole numbers"),
        ({"trigger_duration": [-2.9, 2.9]}, "trigger_duration must be whole"),
        ({"n_components": 2.5}, "n_components must be whole numbers"),
        ({"figure_components": ["x"]}, "figure_components must be whole"),
    ],
)
def test_guards_raise_before_fitting(fluo, behavior, defaults, over, match):
    with pytest.raises(ValueError, match=match):
        _prepare(fluo, behavior, defaults, **over)


def test_string_params_from_hand_edited_yaml_are_coerced(fluo, behavior, defaults):
    prepared = _prepare(
        fluo,
        behavior,
        defaults,
        figure_components=["0", 1],
        feature_columns="0, 2",
        trigger_duration=[-10.0, "10"],
    )
    assert prepared["figure_components"] == [0, 1]
    assert prepared["duration"] == [-10, 10]
    assert prepared["mean"].shape == (N_ROI, 20, 2, 2)


def test_n_components_is_checked_after_iscell_filtering(behavior, defaults, fluo):
    iscell = np.zeros(N_ROI, dtype=int)
    iscell[:5] = 1
    with pytest.raises(ValueError, match=r"n_components \(8\) exceeds .* \(5\)"):
        mod.prepare_inputs(fluo, behavior, iscell, defaults)


class FakeDPCA:
    def __init__(self, labels, regularizer, n_components, n_iter):
        self.n_components = n_components
        self.keys = [
            "".join(c)
            for r in range(1, len(labels) + 1)
            for c in itertools.combinations(labels, r)
        ]

    def fit_transform(self, X):
        self.explained_variance_ratio_ = {
            k: [0.5] * self.n_components for k in self.keys
        }
        return {
            k: np.full((self.n_components,) + X.shape[1:], float(i))
            for i, k in enumerate(self.keys)
        }


@pytest.fixture
def fake_dpca(monkeypatch):
    module = types.ModuleType("dPCA")
    module.dPCA = types.SimpleNamespace(dPCA=FakeDPCA)
    monkeypatch.setitem(sys.modules, "dPCA", module)


def test_dpca_fit_assembles_outputs_with_stub_library(
    fluo, behavior, defaults, fake_dpca
):
    rng_before = np.random.get_state()[1].copy()
    info = mod.dpca_fit(
        FluoData(fluo), BehaviorData(behavior), OUTPUT_DIR, params=defaults
    )
    assert np.array_equal(np.random.get_state()[1], rng_before)

    assert [k for k in info if k != "nwbfile"] == HEATMAP_NAMES
    heat = info["c-component1"]
    assert heat.data.shape == (4, 20)
    assert list(heat.columns) == list(range(-10, 10))
    assert np.all(heat.data == 2.0)
    post = info["nwbfile"][NWBDATASET.POSTPROCESS]["dpca_x"]
    assert post["t"].shape == (8, 20, 2, 2)
    assert post["explained_variance_ratio_tbc"].shape == (8,)
    assert post["levels_b"].tolist() == [0.0, 1.0]
    assert post["levels_c"].tolist() == [0.0, 1.0]
    assert post["trial_counts"].sum() == 99


@pytest.mark.skipif(
    importlib.util.find_spec("dPCA") is None, reason="dPCA only in the conda env"
)
def test_dpca_fit_end_to_end_is_reproducible(fluo, behavior, defaults):
    def run():
        return mod.dpca_fit(
            FluoData(fluo), BehaviorData(behavior), OUTPUT_DIR, params=defaults
        )

    first, second = run(), run()
    assert [k for k in first if k != "nwbfile"] == HEATMAP_NAMES
    assert first["t-component0"].data.shape == (4, 20)
    for name in HEATMAP_NAMES:
        assert np.array_equal(first[name].data, second[name].data)
    post = first["nwbfile"][NWBDATASET.POSTPROCESS]["dpca_x"]
    assert post["explained_variance_ratio_t"].shape == (8,)
    assert post["levels_c"].tolist() == [0.0, 1.0]
