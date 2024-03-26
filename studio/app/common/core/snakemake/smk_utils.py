import os

from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.utils.filepath_finder import find_condaenv_filepath
from studio.app.const import FILETYPE
from studio.app.dir_path import DIRPATH
from studio.app.wrappers import wrapper_dict


class SmkUtils:
    @classmethod
    def input(cls, details):
        if details["type"] in [FILETYPE.IMAGE]:
            return [join_filepath([DIRPATH.INPUT_DIR, x]) for x in details["input"]]
        elif details["type"] in [
            FILETYPE.CSV,
            FILETYPE.BEHAVIOR,
            FILETYPE.HDF5,
            FILETYPE.MATLAB,
            FILETYPE.MICROSCOPE,
        ]:
            return join_filepath([DIRPATH.INPUT_DIR, details["input"]])
        else:
            return [join_filepath([DIRPATH.OUTPUT_DIR, x]) for x in details["input"]]

    @classmethod
    def output(cls, details):
        return join_filepath([DIRPATH.OUTPUT_DIR, details["output"]])

    @classmethod
    def conda(cls, details):
        if details["type"] in [
            FILETYPE.IMAGE,
            FILETYPE.CSV,
            FILETYPE.BEHAVIOR,
            FILETYPE.HDF5,
            FILETYPE.MATLAB,
            FILETYPE.MICROSCOPE,
        ]:
            return None

        wrapper = cls.dict2leaf(wrapper_dict, details["path"].split("/"))

        if "conda_name" in wrapper:
            conda_name = wrapper["conda_name"]
            conda_filepath = f"{DIRPATH.CONDAENV_DIR}/envs/{conda_name}"
            if os.path.exists(conda_filepath):
                return conda_filepath
            else:
                return find_condaenv_filepath(conda_name)

        return None

    @classmethod
    def dict2leaf(cls, root_dict: dict, path_list):
        path = path_list.pop(0)
        if len(path_list) > 0:
            return cls.dict2leaf(root_dict[path], path_list)
        else:
            return root_dict[path]
