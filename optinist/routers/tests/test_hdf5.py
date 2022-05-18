import pytest
import shutil
from fastapi.testclient import TestClient

from optinist.routers.hdf5 import router
from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import create_directory, join_filepath

client = TestClient(router)


def test_hdf5():
    filepath = join_filepath([
        "test_data",
        "test.nwb"
    ])
    output_dir = join_filepath([DIRPATH.INPUT_DIR, "test_data"])
    create_directory(output_dir)
    shutil.copyfile(
        join_filepath([DIRPATH.ROOT_DIR, filepath]),
        join_filepath([output_dir, "test.nwb"])
    )
    response = client.get(f"/hdf5/{filepath}")
    data = response.json()
    assert response.status_code == 200
    assert isinstance(data, list)
    assert isinstance(data[0], dict)

    shutil.rmtree(output_dir)
