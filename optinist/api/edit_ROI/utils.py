import os
from typing import Tuple

from optinist.api.dataclass.base import BaseData
from optinist.api.nwb.nwb_creater import overwrite_nwb


def save_edit_ROI_data(func):
    def wrapper(*args, **kwargs):
        info = func(*args, **kwargs)

        node_dirpath = args[0]
        node_name = os.path.basename(node_dirpath)[:-11]

        for k, v in info.items():
            if isinstance(v, BaseData):
                v.save_json(node_dirpath)

            if k == 'nwbfile':
                overwrite_nwb(v, node_dirpath, f'{node_name}.nwb')

    return wrapper


def create_mask(x: int, y: int, width: int, height: int, dims: Tuple[int, int]):
    import numpy as np

    x, y, width, height = round(x), round(y), round(width), round(height)

    x_coords = np.arange(0, dims[0])
    y_coords = np.arange(0, dims[1])
    xx, yy = np.meshgrid(x_coords, y_coords)

    # Calculate the distance of each pixel from the center of the ellipse
    a = width / 2
    b = height / 2
    distance = ((xx - x) / a) ** 2 + ((yy - y) / b) ** 2

    # Set the pixels within the ellipse to 1, and the pixels outside the ellipse to 0
    ellipse = np.zeros(dims)
    ellipse[distance <= 1] = 1

    return ellipse
