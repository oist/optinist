from optinist.api.pickle.pickle_reader import PickleReader

filepath = "/tmp/optinist/output/0123/func1/func1.pkl"


def test_PickleReader():
    data = PickleReader.read(filepath)

    assert data == "abc"
