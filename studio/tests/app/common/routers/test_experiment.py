import os
import shutil

from fastapi.testclient import TestClient

from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.routers.experiment import router
from studio.app.dir_path import DIRPATH

client = TestClient(router)

unique_id = "0123"

output_test_dir = f"{DIRPATH.DATA_DIR}/output_test"

shutil.copytree(
    f"{output_test_dir}/{unique_id}",
    f"{DIRPATH.OUTPUT_DIR}/{unique_id}",
    dirs_exist_ok=True,
)


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
    dirpath = join_filepath([f"{DIRPATH.DATA_DIR}/output", dirname])
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

    response = client.post("/experiments/delete", json={"uidList": uidList})
    assert response.status_code == 200

    for name in uidList:
        dirpath = join_filepath([DIRPATH.OUTPUT_DIR, name])
        assert not os.path.exists(dirpath)


def test_download_config():
    response = client.get(f"/experiments/download/config/{unique_id}")
    assert response.status_code == 200
    assert response.url == "http://testserver/experiments/download/config/0123"


def test_expt_rename():
    origin_name = "New flow"
    new_name = "TEST RENAME WORKFLOW"
    response = client.patch(
        f"/experiments/{unique_id}/rename", json={"new_name": new_name}
    )
    data = response.json()
    assert response.status_code == 200
    assert data["name"] == new_name
    response = client.patch(
        f"/experiments/{unique_id}/rename", json={"new_name": origin_name}
    )
