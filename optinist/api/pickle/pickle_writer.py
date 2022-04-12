import os
import pickle

from optinist.api.utils.filepath_creater import join_filepath


class PickleWriter:
    @classmethod
    def write(cls, pickle_path, info):
        # ファイル保存先
        os.makedirs(
            join_filepath(pickle_path.split("/")[:-1]),
            exist_ok=True
        )
        with open(pickle_path, 'wb') as f:
            pickle.dump(info, f)
