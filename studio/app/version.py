import os
import re

from studio.app.dir_path import DIRPATH


def get_app_version_from_pyproject() -> str:
    """Extract version from pyproject.toml."""
    # Look for pyproject.toml in parent directories
    pyproject_path = os.path.join(DIRPATH.ROOT_DIR, "pyproject.toml")

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
    APP_VERSION = get_app_version_from_pyproject()
