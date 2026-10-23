import os
import pickle
import traceback
from glob import glob

from studio.app.common.core.utils.filepath_creater import (
    create_directory,
    join_filepath,
)


class PickleReader:
    @classmethod
    def read(cls, filepath):
        with open(filepath, "rb") as f:
            return pickle.load(f)

    @classmethod
    def search_node_pickle_path(cls, search_path: str) -> str:
        paths = list(
            set(glob(join_filepath([search_path, "*.pkl"])))
            - set(glob(join_filepath([search_path, "tmp_*.pkl"])))
        )
        path = paths[0] if paths else None

        return path

    @staticmethod
    def check_is_valid_node_pickle(data):
        """
        Checks whether the node processing result pickle is valid
        - If the processing is successful, the result information is stored in a dict
        - If the processing fails, the error information is stored in a list[str]
        """
        is_valid = (data is not None) and (type(data) is dict)

        return is_valid

    @staticmethod
    def check_is_error_node_pickle(data):
        """
        @see Note on `check_is_valid_node_pickle`
        """
        is_error = (data is None) or isinstance(data, (list, str))

        return is_error


class PickleWriter:
    @classmethod
    def write(cls, pickle_path, info):
        dirpath = join_filepath(pickle_path.split("/")[:-1])
        create_directory(dirpath)

        # Note: Use temporary files to avoid read during file writing.
        tmp_pickle_path = f"{pickle_path}.tmp"

        if os.path.exists(pickle_path):
            os.remove(pickle_path)

        with open(tmp_pickle_path, "wb") as f:
            pickle.dump(info, f)
            f.flush()

        os.rename(tmp_pickle_path, pickle_path)

    @classmethod
    def write_error(cls, pickle_path, err: Exception):
        err_msg = list(traceback.TracebackException.from_exception(err).format())

        cls.write(pickle_path, err_msg)

    @classmethod
    def overwrite(cls, pickle_path, info):
        with open(pickle_path, "rb") as f:
            old_pkl = pickle.load(f)

        if isinstance(old_pkl, dict) and isinstance(info, dict):
            old_pkl.update(info)

            # Note: Use temporary files to avoid read during file writing.
            tmp_pickle_path = f"{pickle_path}.tmp"

            if os.path.exists(pickle_path):
                os.remove(pickle_path)

            with open(tmp_pickle_path, "wb") as f:
                pickle.dump(old_pkl, f)
                f.flush()

            os.rename(tmp_pickle_path, pickle_path)
