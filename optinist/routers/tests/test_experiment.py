import pytest
from fastapi.testclient import TestClient

from optinist.routers.experiment import router

client = TestClient(router)


def test_read_experiment():
    pass