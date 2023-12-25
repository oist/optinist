from typing import Tuple

import numpy as np

from studio.app.optinist.schemas.roi import RoiPos


def create_ellipse_mask(shape: Tuple[int, int], roi_pos: RoiPos):
    x, y, width, height = (
        round(roi_pos.posx),
        round(roi_pos.posy),
        round(roi_pos.sizex),
        round(roi_pos.sizey),
    )

    x_coords = np.arange(0, shape[0])
    y_coords = np.arange(0, shape[1])
    xx, yy = np.meshgrid(x_coords, y_coords)

    # Calculate the distance of each pixel from the center of the ellipse
    a = width / 2
    b = height / 2
    distance = ((xx - x) / a) ** 2 + ((yy - y) / b) ** 2

    # Set the pixels within the ellipse to 1 and the pixels outside to NaN
    ellipse = np.empty(shape)
    ellipse[:] = np.nan
    ellipse[distance <= 1] = 1

    return ellipse
