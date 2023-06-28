import os

from optinist.api.dir_path import DIRPATH
from optinist.api.pickle.pickle_writer import PickleWriter

filepath = f"{DIRPATH.OPTINIST_DIR}/output/0123/func2/func2.pkl"


def test_PickleWriter():
    PickleWriter.write(
        filepath,
        {"abc": None},
    )

    assert os.path.exists(filepath)


def test_PickleWriter_overwrite():
    PickleWriter.overwrite(
        filepath,
        {"test": 1},
    )

    assert os.path.exists(filepath)
