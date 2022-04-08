import os


def join_filepath(path_list):
    if isinstance(path_list, str):
        return path_list
    elif isinstance(path_list, list):
        return "/".join(path_list)
    assert False, "Path is not list"


def create_filepath(dirname, filename):
    if not os.path.exists(dirname):
        os.makedirs(dirname, exist_ok=True)

    return join_filepath([dirname, filename])
