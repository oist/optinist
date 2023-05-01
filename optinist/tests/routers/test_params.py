import pytest
from fastapi.testclient import TestClient

from optinist.routers.params import router

client = TestClient(router)


def test_params():
    response = client.get("/params/caiman_mc")
    data = response.json()

    assert response.status_code == 200
    assert isinstance(data, dict)

    assert isinstance(data["border_nan"], str)
    assert data["border_nan"] == "copy"

    assert isinstance(data["use_cuda"], bool)
    assert data["use_cuda"] == False


def test_snakemake_params():
    response = client.get("/snakemake")
    data = response.json()

    assert response.status_code == 200
    assert isinstance(data, dict)

    assert isinstance(data["cores"], int)
    assert data["cores"] == 2

    assert isinstance(data["use_conda"], bool)
    assert data["use_conda"] == True


def test_nwb_params():
    response = client.get("/nwb")
    data = response.json()

    assert response.status_code == 200
    assert isinstance(data, dict)

    assert isinstance(data["session_description"], str)
    assert data["session_description"] == "optinist"
