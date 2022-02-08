import os
import copy
import numpy as np


def get_file_path(_dir, file_name):
    if not os.path.exists(_dir):
        os.makedirs(_dir, exist_ok=True)

    return join_file_path([_dir, f'{file_name}.json'])


def get_images_list(data):
    if len(data.shape) == 2:
        save_data = copy.deepcopy(data)
    elif len(data.shape) == 3:
        save_data = copy.deepcopy(data[:10])
    else:
        assert False, 'data is error'

    if len(data.shape) == 2:
        data = data[np.newaxis, :, :]
        save_data = save_data[np.newaxis, :, :]

    images = []
    for i, _img in enumerate(save_data):
        images.append(_img.tolist())
    
    return images


def join_file_path(path_list):
    if isinstance(path_list, str):
        return path_list
    elif isinstance(path_list, list):
        return "/".join(path_list)
    assert False, "Path is not list"
