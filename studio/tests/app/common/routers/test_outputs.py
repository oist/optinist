from studio.app.dir_path import DIRPATH

timeseries_dirpath = f"{DIRPATH.DATA_DIR}/output/default/0123/func1/fluorescence.json"


def test_inittimedata(client):
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


def test_timedata(client):
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


def test_alltimedata(client):
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
workspace_id = "1"


def test_image(client):
    response = client.get(f"/outputs/image/{tif_filepath}?workspace_id={workspace_id}")
    data = response.json()

    assert response.status_code == 200
    assert isinstance(data, dict)
