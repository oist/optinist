import os
import pandas as pd
import cv2
import numpy as np
from PIL import Image, ImageSequence


def save_tiff_to_json(tiff_file_path, maxidx=10):
    print("save_tiff_to_json")
    folder_path = os.path.dirname(tiff_file_path)
    file_name, ext = os.path.splitext(os.path.basename(tiff_file_path))

    # Tiff画像を読み込む
    tiffs = []
    image = Image.open(tiff_file_path)

    for i, page in enumerate(ImageSequence.Iterator(image)):
        # print(page)
        tiffs.append(np.array(page.resize((150, 150))))

        if i >= maxidx:
            break

    tiffs = np.array(tiffs)

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
