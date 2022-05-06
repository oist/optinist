import pytest
from fastapi.testclient import TestClient
import shutil
import os

from optinist.routers.experiment import router
from optinist.api.utils.filepath_creater import create_directory, join_filepath
from optinist.api.dir_path import DIRPATH

client = TestClient(router)


def test_read_experiment():
    output_dir = join_filepath([DIRPATH.OUTPUT_DIR, "test"]) 
    create_directory(output_dir)
    shutil.copy(
        join_filepath([DIRPATH.ROOT_DIR, "test_data", "experiment.yaml"]),
        join_filepath([output_dir, "experiment.yaml"])
    )
    response = client.get("/experiments")
    data = response.json()
    assert response.status_code == 200
    assert isinstance(data, dict)
    assert isinstance(data[next(iter(data))], dict)

    shutil.rmtree(output_dir)


def test_delete_experiment():
    dirpath = join_filepath([DIRPATH.OUTPUT_DIR, "aa"])
    os.makedirs(dirpath, exist_ok=True)
    assert os.path.exists(dirpath)
    response = client.delete("/experiments/aa")
    assert response.status_code == 200
    assert not os.path.exists(dirpath)


def test_delete_experiment_list():
    uidList = ["aa", "bb"]
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


def test_download_experiment():
    dirpath = join_filepath([DIRPATH.OPTINIST_DIR, "aa", "bb"])
    os.makedirs(dirpath, exist_ok=True)
    with open(join_filepath([dirpath, "aa.nwb"]), "w") as f:
        f.write("aa")

    response = client.get("/experiments/download/nwb/aa")

    assert response.status_code == 200
