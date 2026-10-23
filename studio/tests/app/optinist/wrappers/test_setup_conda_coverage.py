import os

import pytest
import yaml

from studio.app.dir_path import DIRPATH
from studio.app.wrappers import wrapper_dict

SETUP_CONDA = wrapper_dict["maintenance"]["setup_conda"]
SAMPLE_OUTPUT_DIR = os.path.join(DIRPATH.ROOT_DIR, "sample_data/maintenance/output")


def leaf_wrappers(tree):
    for value in tree.values():
        if "function" in value:
            yield value
        else:
            yield from leaf_wrappers(value)


CONDA_NAMES = sorted(
    {w["conda_name"] for w in leaf_wrappers(wrapper_dict) if w.get("conda_name")}
)


@pytest.mark.parametrize("conda_name", CONDA_NAMES)
def test_every_conda_env_can_be_created_from_the_dialog(conda_name):
    # the "Create env" dialog reproduces and runs the sample named setup_conda_<env>
    uid = f"setup_conda_{conda_name}"
    assert SETUP_CONDA[uid]["conda_name"] == conda_name

    with open(os.path.join(SAMPLE_OUTPUT_DIR, uid, "experiment.yaml")) as f:
        assert yaml.safe_load(f)["unique_id"] == uid
