import pytest
import os

from optinist.cui_api.dir_path import DIRPATH

def test_dir_path():
    assert os.path.exists(DIRPATH.BASE_DIR)
    assert os.path.exists(DIRPATH.ROOT_DIR)
