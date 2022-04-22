import os

from optinist.api.dir_path import DIRPATH


def join_filepath(path_list):
    if isinstance(path_list, str):
        return path_list
    elif isinstance(path_list, list):
        return "/".join(path_list)
    assert False, "Path is not list"


def create_filepath(dirname, filename):
    create_directory(dirname)

    return join_filepath([dirname, filename])


def get_pickle_file(unique_id, node_id, algo_name):
    return join_filepath([
        DIRPATH.BASE_DIR,
        unique_id,
        node_id,
        f"{algo_name}.pkl"
    ])


def create_directory(dirpath):
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
