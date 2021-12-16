import os
import imageio
import pandas as pd
import cv2


def save_tiff_to_json(tiff_file_path, maxidx=10):
    folder_path = os.path.dirname(tiff_file_path)
    file_name, ext = os.path.splitext(os.path.basename(tiff_file_path))

    tiffs = imageio.volread(tiff_file_path)[:maxidx]
    # import pdb; pdb.set_trace()
    tiffs = cv2.resize(tiffs.transpose(1, 2, 0), (150, 150)).transpose(2, 0, 1)

    images = []
    for i, _img in enumerate(tiffs):
        images.append(_img.tolist())

    pd.DataFrame(images).to_json(
        os.path.join(folder_path, f'{file_name}.json'), indent=4, orient="values")

def save_csv_to_json(csv_file_path):
    folder_path = os.path.dirname(csv_file_path)
    file_name, ext = os.path.splitext(os.path.basename(csv_file_path))
    pd.read_csv(csv_file_path).to_json(
        os.path.join(folder_path, f'{file_name}.json'), indent=4,orient="split")
