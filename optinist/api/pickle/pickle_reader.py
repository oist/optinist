import pickle


class PickleReader:
    @classmethod
    def read(cls, filepath):
        with open(filepath, "rb") as f:
            return pickle.load(f)
