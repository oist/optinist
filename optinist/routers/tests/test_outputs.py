from fastapi.testclient import TestClient

from optinist.routers.outputs import router
from optinist.cui_api.dir_path import DIRPATH
from optinist.cui_api.filepath_creater import join_filepath

client = TestClient(router)


def test_read_file():
    dirpath = join_filepath([DIRPATH.ROOT_DIR, "test_data", "fluorescence.json"])
    response = client.get(f"/outputs/timedata/{dirpath}/?index=0")
    data = response.json()
    print(data)
    assert response.status_code == 200
    assert isinstance(data, dict)
    assert isinstance(data["data"], dict)
    assert isinstance(data["std"], dict)
    assert isinstance(data["xrange"], list)
    assert data["data"]["0"]["0"] == 479.916595459
