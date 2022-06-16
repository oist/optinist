import pytest
from fastapi.testclient import TestClient
import os

from optinist.routers.experiment import router
from optinist.api.utils.filepath_creater import join_filepath
from optinist.api.dir_path import DIRPATH

client = TestClient(router)

unique_id = "0123"


def test_get():
    response = client.get("/experiments")
    data = response.json()

    assert response.status_code == 200
    assert isinstance(data, dict)
    assert isinstance(data[next(iter(data))], dict)


def test_import():
    response = client.get(f"/experiments/import/{unique_id}")
    data = response.json()

    assert response.status_code == 200
    assert isinstance(data, dict)


def test_delete():
    dirname = "delete_dir"
    dirpath = join_filepath([DIRPATH.OUTPUT_DIR, dirname])
    os.makedirs(dirpath, exist_ok=True)
    assert os.path.exists(dirpath)
    response = client.delete(f"/experiments/{dirname}")
    assert response.status_code == 200
    assert not os.path.exists(dirpath)


def test_delete_list():
    uidList = ["delete_dir1", "delete_dir2"]
    for name in uidList:
        dirpath = join_filepath([DIRPATH.OUTPUT_DIR, name])
        os.makedirs(dirpath, exist_ok=True)
        assert os.path.exists(dirpath)

    response = client.post(
        "/experiments/delete", json={"uidList": uidList}
    )
    assert response.status_code == 200

    for name in uidList:
        dirpath = join_filepath([DIRPATH.OUTPUT_DIR, name])
        assert not os.path.exists(dirpath)


def test_download_nwb():
    response = client.get(f"/experiments/download/nwb/{unique_id}")

    assert response.status_code == 200
    assert response.url == "http://testserver/experiments/download/nwb/0123"


def test_download_nwb_function():
    function_id = "func1"
    response = client.get(f"/experiments/download/nwb/{unique_id}/{function_id}")
    assert response.status_code == 200
    assert response.url == "http://testserver/experiments/download/nwb/0123/func1"


def test_download_config():
    response = client.get(f"/experiments/download/config/{unique_id}")
    assert response.status_code == 200
    assert response.url == "http://testserver/experiments/download/config/0123"
