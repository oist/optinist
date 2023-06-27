import os

from studio.core.dir_path import DIRPATH
from studio.core.utils.pickle_handler import PickleReader, PickleWriter

filepath = f"{DIRPATH.OPTINIST_DIR}/output/0123/func2/func2.pkl"


def test_PickleWriter():
    PickleWriter.write(
        filepath,
        "abc",
    )

    assert os.path.exists(filepath)


filepath = f"{DIRPATH.OPTINIST_DIR}/output_test/0123/func1/func1.pkl"


def test_PickleReader():
    data = PickleReader.read(filepath)

    assert data == "abc"
