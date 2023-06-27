import os

from optinist.services.dir_path import DIRPATH


def test_dir_path():
    assert os.path.exists(DIRPATH.OPTINIST_DIR)
    assert os.path.exists(DIRPATH.ROOT_DIR)
