from fastapi.testclient import TestClient

from studio.app.optinist.routers.nwb import router

client = TestClient(router)


def test_nwb_params():
    response = client.get("/nwb")
    data = response.json()

    assert response.status_code == 200
    assert isinstance(data, dict)

    assert isinstance(data["session_description"], str)
    assert data["session_description"] == "optinist"
