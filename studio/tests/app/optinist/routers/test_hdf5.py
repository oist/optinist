from studio.app.dir_path import DIRPATH
from studio.app.optinist.routers.hdf5 import HDF5Getter
from studio.app.optinist.schemas.hdf5 import HDF5Node

input_filepath = "files/test.nwb"
workspace_id = "1"


def test_hdf5(client):
    response = client.get(f"/hdf5/{input_filepath}?workspace_id={workspace_id}")
    data = response.json()

    assert response.status_code == 200
    assert isinstance(data, list)
    assert isinstance(data[0], dict)


def test_HDF5Getter():
    output = HDF5Getter.get(f"{DIRPATH.DATA_DIR}/input/1/files/test.nwb")

    assert isinstance(output, list)
    assert isinstance(output[0], HDF5Node)
