import os
import re
from importlib.metadata import version as package_version

from studio.app.dir_path import DIRPATH

PYPROJECT_PATH = os.path.join(DIRPATH.ROOT_DIR, "pyproject.toml")


def get_app_version() -> str:
    """Get version from pyproject.toml in a source tree, else from package metadata."""
    # pip installs ship no pyproject.toml; editable installs keep stale metadata
    if os.path.exists(PYPROJECT_PATH):
        return get_app_version_from_pyproject(PYPROJECT_PATH)
    return package_version("optinist")


def get_app_version_from_pyproject(pyproject_path: str) -> str:
    """Extract version from pyproject.toml."""
    # Read and parse the version from pyproject.toml
    with open(pyproject_path, "r") as f:
        content = f.read()

    # Use regex to extract version
    version_match = re.search(r'version\s*=\s*["\']([^"\']+)["\']', content)
    if version_match:
        version = version_match.group(1)
        return version
    else:
        return "1.0.0"


class Version:
    APP_VERSION = get_app_version()
