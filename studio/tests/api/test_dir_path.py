import pytest
import os

from studio.api.dir_path import DIRPATH

def test_dir_path():
    assert os.path.exists(DIRPATH.studio_DIR)
    assert os.path.exists(DIRPATH.ROOT_DIR)
