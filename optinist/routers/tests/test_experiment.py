import pytest
from fastapi.testclient import TestClient


client = TestClient(router)


def test_read_experiment():
    pass