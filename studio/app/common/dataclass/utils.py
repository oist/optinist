import copy

import numpy as np


def create_images_list(data):
    assert len(data.shape) == 2, "data is error"

    save_data = copy.deepcopy(data)
    data = data[np.newaxis, :, :]
    save_data = save_data[np.newaxis, :, :]

    images = []
    for _img in save_data:
        images.append(_img.tolist())

    return images


def save_thumbnail(plot_file):
    # Note: In the barebone-studio version, it does nothing.
    pass
