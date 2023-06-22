import pickle

from optinist.api.utils.filepath_creater import create_directory, join_filepath


class PickleWriter:
    @classmethod
    def write(cls, pickle_path, info):
        # ファイル保存先
        dirpath = join_filepath(pickle_path.split("/")[:-1])
        create_directory(dirpath)
        with open(pickle_path, "wb") as f:
            pickle.dump(info, f)

    @classmethod
    def overwrite(cls, pickle_path, info):
        with open(pickle_path, "rb") as f:
            old_pkl = pickle.load(f)

        if type(old_pkl) == type(info) == dict:
            old_pkl.update(info)

            with open(pickle_path, "wb") as f:
                pickle.dump(old_pkl, f)
