import os

from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import join_filepath
from optinist.routers.model import FILETYPE
from optinist.wrappers import wrapper_dict


class SmkUtils:
    @classmethod
    def input(cls, details):
        if details["type"] in [FILETYPE.IMAGE]:
            return [join_filepath([DIRPATH.INPUT_DIR, x]) for x in details["input"]]
        elif details["type"] in [FILETYPE.CSV, FILETYPE.BEHAVIOR, FILETYPE.HDF5]:
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
        ]:
            return None

        wrapper = cls.dict2leaf(wrapper_dict, details["path"].split("/"))

        if "conda_name" in wrapper:
            _filename = wrapper["conda_name"]
            conda_filepath = f"{DIRPATH.CONDAENV_DIR}/envs/{_filename}"
            if os.path.exists(conda_filepath):
                return conda_filepath

        if "conda_yaml" in wrapper:
            conda_yaml = wrapper["conda_yaml"]
            return f"{DIRPATH.CONDAYML_DIR}/{conda_yaml}"

        return None

    @classmethod
    def dict2leaf(cls, root_dict: dict, path_list):
        path = path_list.pop(0)
        if len(path_list) > 0:
            return cls.dict2leaf(root_dict[path], path_list)
        else:
            return root_dict[path]
