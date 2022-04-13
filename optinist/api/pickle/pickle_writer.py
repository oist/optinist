import os
import pickle

from optinist.api.utils.filepath_creater import join_filepath


class PickleWriter:
    @classmethod
    def write(cls, pickle_path, info):
        # ファイル保存先
        dirpath = join_filepath(pickle_path.split("/")[:-1])
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)
        with open(pickle_path, 'wb') as f:
            pickle.dump(info, f)
