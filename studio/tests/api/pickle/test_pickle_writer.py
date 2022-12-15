import pytest
import os

from studio.api.pickle.pickle_writer import PickleWriter


filepath = "/tmp/studio/output/0123/func2/func2.pkl"


def test_PickleWriter():
    PickleWriter.write(
        filepath,
        "abc",
    )

    assert os.path.exists(filepath)
