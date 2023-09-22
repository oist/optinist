from studio.app.dir_path import DIRPATH

workspace_id = "default"
unique_id = "0123"
output_test_dir = f"{DIRPATH.DATA_DIR}/output_test"


def test_nwb_params(client):
    response = client.get("/nwb")
    data = response.json()

    assert response.status_code == 200
    assert isinstance(data, dict)

    assert isinstance(data["session_description"], str)
    assert data["session_description"] == "optinist"


def test_download_nwb(client):
    response = client.get(f"/experiments/download/nwb/{workspace_id}/{unique_id}")

    assert response.status_code == 200
    assert response.url == "http://testserver/experiments/download/nwb/default/0123"


def test_download_nwb_function(client):
    function_id = "func1"
    response = client.get(
        f"/experiments/download/nwb/{workspace_id}/{unique_id}/{function_id}"
    )
    assert response.status_code == 200
    assert (
        response.url == "http://testserver/experiments/download/nwb/default/0123/func1"
    )
