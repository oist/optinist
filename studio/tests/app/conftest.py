import shutil
from typing import Generator

import pytest
from fastapi.testclient import TestClient

from studio.__main_unit__ import app, skip_dependencies
from studio.app.common.core.auth.auth_dependencies import (
    get_admin_user,
    get_current_user,
)
from studio.app.common.core.workspace.workspace_dependencies import (
    is_workspace_available,
    is_workspace_owner,
)
from studio.app.dir_path import DIRPATH


@pytest.fixture(scope="session", autouse=True)
def session_fixture():
    app.dependency_overrides[get_admin_user] = skip_dependencies
    app.dependency_overrides[get_current_user] = skip_dependencies
    app.dependency_overrides[is_workspace_available] = skip_dependencies
    app.dependency_overrides[is_workspace_owner] = skip_dependencies

    yield

    shutil.rmtree(f"{DIRPATH.DATA_DIR}/output")


@pytest.fixture(scope="module")
def client() -> Generator:
    with TestClient(app) as c:
        yield c
