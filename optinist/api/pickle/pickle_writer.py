import pickle

from optinist.api.utils.filepath_creater import create_directory, join_filepath


class PickleWriter:
    @classmethod
    def write(cls, pickle_path, info):
        # ファイル保存先
        dirpath = join_filepath(pickle_path.split("/")[:-1])
        create_directory(dirpath)
        with open(pickle_path, 'wb') as f:
            pickle.dump(info, f)
