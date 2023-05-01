import pytest
from fastapi.testclient import TestClient

from optinist.routers.hdf5 import HDF5Getter, router
from optinist.routers.model import HDF5Node

client = TestClient(router)

input_filepath = "files/test.nwb"


def test_hdf5():
    response = client.get(f"/hdf5/{input_filepath}")
    data = response.json()

    assert response.status_code == 200
    assert isinstance(data, list)
    assert isinstance(data[0], dict)


def test_HDF5Getter():
    output = HDF5Getter.get("/tmp/optinist/input/files/test.nwb")

    assert isinstance(output, list)
    assert isinstance(output[0], HDF5Node)
