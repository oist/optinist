import pytest
import os

from optinist.cui_api.dir_path import DIRPATH
from optinist.cui_api.filepath_creater import (
    create_filepath,
    join_filepath
)

def test_create_filepath():
    filepath = create_filepath(
        join_filepath([DIRPATH.ROOT_DIR, "test_data"]),
        "json_filepath.json",
    )

    with open(filepath, "w") as f:
        f.write("")

    assert isinstance(filepath, str)
    assert os.path.exists(filepath)
