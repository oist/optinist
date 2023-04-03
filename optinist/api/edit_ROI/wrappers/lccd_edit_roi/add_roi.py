import os
from typing import Tuple

from optinist.api.dataclass.dataclass import *


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


def execute_add_ROI(node_dirpath, posx, posy, sizex, sizey):
    import numpy as np


    lccd_data = np.load(os.path.join(node_dirpath, 'lccd.npy'), allow_pickle=True).item()

    images = lccd_data.get('images', None)
    roi = lccd_data.get('roi', None)
    shape = images.shape[:2]
    num_frames = images.shape[2]
    dff_f0_frames = 100
    dff_f0_percentile = 8

    new_roi = create_mask(posx, posy, sizex, sizey, shape).reshape(-1, 1)
    roi = np.hstack((new_roi,roi))


    num_cell = roi.shape[1]
    images = images.reshape([images.shape[0] * images.shape[1], images.shape[2]])
    timeseries = np.zeros([num_cell, num_frames])

    # Get ROI list and extract Fluor
    roi_list = []
    for i in range(num_cell):
        roi_list.append((roi[:, i].reshape(shape)) * (i + 1))
        timeseries[i, :] = np.mean(images[roi[:, i] > 0, :], axis=0)
    roi_list = np.stack(roi_list)
    roi_list[roi_list == 0] = np.nan

    # Get DFF
    timeseries_dff = np.ones([num_cell, num_frames]) * np.nan
    for i in range(num_cell):
        for k in range(num_frames):
            if (k - dff_f0_frames >= 0) and (k + dff_f0_frames < num_frames):
                f0 = np.percentile(
                    timeseries[i, k - dff_f0_frames : k + dff_f0_frames],
                    dff_f0_percentile,
                )
                timeseries_dff[i, k] = (timeseries[i, k] - f0) / f0

    lccd_data['roi'] = roi

    info = {
        'lccd': LccdData(lccd_data),
        'cell_roi': RoiData(np.nanmax(roi_list, axis=0), file_name='cell_roi'),
        'fluorescence': FluoData(timeseries, file_name='fluorescence'),
        'dff': FluoData(timeseries_dff, file_name='dff'),
    }

    for v in info.values():
        if isinstance(v, BaseData):
            v.save_json(node_dirpath)
