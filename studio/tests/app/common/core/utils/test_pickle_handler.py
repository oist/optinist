import os

from studio.app.common.core.utils.pickle_handler import PickleReader, PickleWriter
from studio.app.dir_path import DIRPATH

workspace_id = "default"
unique_id = "pickle_test"

filepath = f"{DIRPATH.OUTPUT_DIR}/{workspace_id}/{unique_id}/func2/func2.pkl"


def test_PickleWriter():
    PickleWriter.write(
        filepath,
        "abc",
    )

    assert os.path.exists(filepath)


filepath = f"{DIRPATH.DATA_DIR}/output_test/{workspace_id}/{unique_id}/func1/func1.pkl"


def test_PickleReader():
    data = PickleReader.read(filepath)

    assert data == "abc"
