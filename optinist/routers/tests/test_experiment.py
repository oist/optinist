import pytest
from fastapi.testclient import TestClient
import os

from optinist.routers.experiment import router, DeleteItem
from optinist.cui_api.experiment_config import ExpConfig
from optinist.cui_api.filepath_creater import join_filepath
from optinist.cui_api.dir_path import DIRPATH

client = TestClient(router)


def test_read_experiment():
    response = client.get("/experiments")
    data = response.json()
    assert response.status_code == 200
    assert isinstance(data, dict)
    assert isinstance(data[next(iter(data))], dict)


def test_delete_experiment():
    dirpath = join_filepath([DIRPATH.BASE_DIR, "aa"])
    os.makedirs(dirpath, exist_ok=True)
    assert os.path.exists(dirpath)
    response = client.delete("/experiments/aa")
    data = response.json()
    assert response.status_code == 200
    assert not os.path.exists(dirpath)


def test_delete_experiment_list():
    uidList = ["aa", "bb"]
    for name in uidList:
        dirpath = join_filepath([DIRPATH.BASE_DIR, name])
        os.makedirs(dirpath, exist_ok=True)
        assert os.path.exists(dirpath)

    response = client.post(
        "/experiments/delete", json={"uidList": uidList}
    )
    assert response.status_code == 200

    for name in uidList:
        dirpath = join_filepath([DIRPATH.BASE_DIR, name])
        assert not os.path.exists(dirpath)


def test_download_experiment():
    dirpath = join_filepath([DIRPATH.BASE_DIR, "aa", "bb"])
    os.makedirs(dirpath, exist_ok=True)
    with open(join_filepath([dirpath, "aa.nwb"]), "w") as f:
        f.write("aa")

    response = client.get("/experiments/download/aa")

    assert response.status_code == 200
