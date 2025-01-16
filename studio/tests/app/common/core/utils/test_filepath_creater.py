import os

from studio.app.dir_path import DIRPATH

workspace_id = "default"

filepath = f"{DIRPATH.OUTPUT_DIR}/{workspace_id}/test.txt"


def test_create_filepath():
    with open(filepath, "w") as f:
        f.write("abc")

    assert isinstance(filepath, str)
    assert os.path.exists(filepath)
