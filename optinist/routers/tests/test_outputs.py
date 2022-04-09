import pytest
import shutil
from fastapi.testclient import TestClient

from optinist.routers.outputs import router

client = TestClient(router)


def test_outputs():
    filepath = join_filepath([DIRPATH.ROOT_DIR, "test_data", "test.nwb"])
    response = client.get(f"/hdf5/{filepath}")
    data = response.json()
    assert response.status_code == 200
    assert isinstance(data, list)
    assert isinstance(data[0], dict)
