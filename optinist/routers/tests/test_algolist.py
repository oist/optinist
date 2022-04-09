import pytest
from fastapi.testclient import TestClient

from optinist.routers.algolist import router


client = TestClient(router)


def test_read_run():
    response = client.get("/algolist")
    data = response.json()
    assert response.status_code == 200
    assert isinstance(data, dict)
    assert "caiman" in data
    assert "children" in data["caiman"]
    assert "caiman_mc" in data["caiman"]["children"]
    assert "args" in data["caiman"]["children"]["caiman_mc"]
    assert "path" in data["caiman"]["children"]["caiman_mc"]
    assert "suite2p" in data
