from fastapi.testclient import TestClient

from studio.app.common.routers.workflow import router

client = TestClient(router)

workspace_id = "default"
unique_id = "0123"


def test_import():
    response = client.get(f"/workflow/import/{workspace_id}/{unique_id}")
    data = response.json()

    assert response.status_code == 200
    assert isinstance(data, dict)
