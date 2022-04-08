import os
import copy
import numpy as np


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