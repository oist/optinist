from fastapi.testclient import TestClient

from studio.app.common.routers.outputs import router
from studio.app.dir_path import DIRPATH

client = TestClient(router)

timeseries_dirpath = f"{DIRPATH.DATA_DIR}/output/default/0123/func1/fluorescence.json"


def test_inittimedata():
    response = client.get(f"/outputs/inittimedata/{timeseries_dirpath}")
    data = response.json()

    assert response.status_code == 200
    assert isinstance(data, dict)
    assert isinstance(data["data"], dict)
    assert isinstance(data["std"], dict)
    assert isinstance(data["xrange"], list)

    assert len(data["data"]) == 67
    assert len(data["data"]["0"]) > 1

    for key, value in data["data"].items():
        if key == "0":
            assert len(value) == 1000
        else:
            assert len(value) == 1


def test_timedata():
    index = 0
    response = client.get(f"/outputs/timedata/{timeseries_dirpath}/?index={index}")
    data = response.json()

    assert response.status_code == 200
    assert isinstance(data, dict)
    assert isinstance(data["data"], dict)
    assert isinstance(data["std"], dict)
    assert isinstance(data["xrange"], list)

    assert str(index) in data["data"]
    assert data["data"]["0"]["0"] == 479.916595459

    index = 1
    response = client.get(f"/outputs/timedata/{timeseries_dirpath}/?index={index}")
    data = response.json()

    assert response.status_code == 200

    assert str(index) in data["data"]
    assert data["data"]["1"]["0"] == 488.6315612793


def test_alltimedata():
    response = client.get(f"/outputs/alltimedata/{timeseries_dirpath}")
    data = response.json()

    assert response.status_code == 200
    assert isinstance(data, dict)
    assert isinstance(data["data"], dict)
    assert isinstance(data["std"], dict)
    assert isinstance(data["xrange"], list)

    assert len(data["data"]) == 67

    for value in data["data"].values():
        assert len(value) == 1000


tif_filepath = "test.tif"


def test_image():
    response = client.get(f"/outputs/image/{tif_filepath}")
    data = response.json()

    assert response.status_code == 200
    assert isinstance(data, dict)
