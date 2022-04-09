import os
import numpy as np
import pandas as pd
from datetime import datetime
from pynwb import NWBHDF5IO
from PIL import Image, ImageSequence
import tifffile
from optinist.cui_api.filepath_creater import join_filepath


def save_tiff2json(tiff_file_path, start_index=None, end_index=None):
    folder_path = os.path.dirname(tiff_file_path)
    file_name, ext = os.path.splitext(os.path.basename(tiff_file_path))

    # Tiff画像を読み込む
    tiffs = []
    image = tifffile.imread(tiff_file_path)
    if image.ndim == 2:
        image = image[np.newaxis, :, :]

    for i, page in enumerate(image):
        if i < start_index-1:
            continue

        if i >= end_index:
            break

        page = np.array(page)
        tiffs.append(page)

    tiffs = np.array(tiffs)

    images = []
    for i, _img in enumerate(tiffs):
        images.append(_img.tolist())

    pd.DataFrame(images).to_json(
        join_filepath([folder_path, f'{file_name}_{str(start_index)}_{str(end_index)}.json']),
        indent=4,
        orient="values"
    )


def save_csv2json(csv_file_path):
    folder_path = os.path.dirname(csv_file_path)
    file_name, ext = os.path.splitext(os.path.basename(csv_file_path))
    pd.read_csv(csv_file_path, header=None).to_json(
        join_filepath([folder_path, f'{file_name}.json']), indent=4, orient='values')
