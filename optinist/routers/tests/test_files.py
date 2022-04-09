import pytest
from fastapi.testclient import TestClient

from optinist.routers.files import router

client = TestClient(router)


def test_create_files():
    response = client.get("/files?file_type=image")
    data = response.json()
    assert response.status_code == 200
    assert isinstance(data, list)
