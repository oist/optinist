from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import join_filepath
from optinist.routers.model import FILETYPE
from optinist.wrappers import wrapper_dict


def smk_func_input(details):
    if details["type"] in [FILETYPE.IMAGE]:
        return [
            join_filepath([DIRPATH.INPUT_DIR, x])
            for x in details["input"]
        ]
    elif details["type"] in [FILETYPE.CSV, FILETYPE.BEHAVIOR, FILETYPE.HDF5]:
        return join_filepath([DIRPATH.INPUT_DIR, details["input"]])
    else:
        return [
            join_filepath([DIRPATH.OUTPUT_DIR, x])
            for x in details["input"]
        ]


def smk_func_output(details):
    return join_filepath([
        DIRPATH.OUTPUT_DIR,
        details["output"]
    ])


def get_conda(details):
    if details["type"] in [FILETYPE.IMAGE, FILETYPE.CSV, FILETYPE.BEHAVIOR, FILETYPE.HDF5]:
        return None

    wrapper = _dict2leaf(
        wrapper_dict,
        details["path"].split('/')
    )

    if "conda" in wrapper:
        env_filename = wrapper['conda']
        return f"{DIRPATH.ROOT_DIR}/rules/envs/{env_filename}"
    else:
        return None


def _dict2leaf(root_dict: dict, path_list):
    path = path_list.pop(0)
    if len(path_list) > 0:
        return _dict2leaf(root_dict[path], path_list)
    else:
        return root_dict[path]
