import pytest
import shutil
from fastapi.testclient import TestClient

from optinist.routers.hdf5 import router
from optinist.cui_api.dir_path import DIRPATH
from optinist.cui_api.filepath_creater import join_filepath

client = TestClient(router)


def test_hdf5():
    filepath = join_filepath([DIRPATH.ROOT_DIR, "test_data", "test.nwb"])
    response = client.get(f"/hdf5/{filepath}")
    data = response.json()
    assert response.status_code == 200
    assert isinstance(data, list)
    assert isinstance(data[0], dict)
