import os

from studio.app.common.core.utils.pickle_handler import PickleReader, PickleWriter
from studio.app.dir_path import DIRPATH

filepath = f"{DIRPATH.DATA_DIR}/output/0123/func2/func2.pkl"


def test_PickleWriter():
    PickleWriter.write(
        filepath,
        "abc",
    )

    assert os.path.exists(filepath)


filepath = f"{DIRPATH.DATA_DIR}/output_test/0123/func1/func1.pkl"


def test_PickleReader():
    data = PickleReader.read(filepath)

    assert data == "abc"
