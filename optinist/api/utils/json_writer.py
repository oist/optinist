import os

import numpy as np
import pandas as pd
import tifffile

from optinist.api.utils.filepath_creater import create_directory, join_filepath


class JsonWriter:
    @classmethod
    def write(cls, filepath, data):
        pd.DataFrame(data).to_json(filepath, indent=4)

    @classmethod
    def write_as_split(cls, filepath, data):
        pd.DataFrame(data).to_json(filepath, indent=4, orient="split")


def save_tiff2json(tiff_filepath, save_dirpath, start_index=None, end_index=None):
    # Tiff画像を読み込む
    tiffs = []
    image = tifffile.imread(tiff_filepath)
    if image.ndim == 2:
        image = image[np.newaxis, :, :]

    for i, page in enumerate(image):
        if i < start_index - 1:
            continue

        if i >= end_index:
            break

        tiffs.append(page.tolist())

    filename, _ = os.path.splitext(os.path.basename(tiff_filepath))
    create_directory(save_dirpath)

    JsonWriter.write_as_split(
        join_filepath(
            [save_dirpath, f"{filename}_{str(start_index)}_{str(end_index)}.json"]
        ),
        pd.DataFrame(tiffs),
    )
