from typing import Tuple


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
