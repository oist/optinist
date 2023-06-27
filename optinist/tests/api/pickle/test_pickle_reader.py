from optinist.services.dir_path import DIRPATH
from optinist.services.pickle.pickle_reader import PickleReader

filepath = f"{DIRPATH.OPTINIST_DIR}/output_test/0123/func1/func1.pkl"


def test_PickleReader():
    data = PickleReader.read(filepath)

    assert data == "abc"
