import shutil

import pytest

from studio.app.dir_path import DIRPATH


@pytest.fixture(scope="session", autouse=True)
def cleanup_testdata():
    yield
    shutil.rmtree(f"{DIRPATH.DATA_DIR}/output")
