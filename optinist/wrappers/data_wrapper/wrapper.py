import os
import numpy as np
import pandas as pd
import imageio
from PIL import Image
import cv2
import tifffile
import copy


class ImageData:
    def __init__(self, data, func_name='image', file_name='image'):
        if type(data) == str:
            self.path = data
            self.data = np.array(imageio.volread(self.path))
            self.json_path = None
        else:
            self.data = data
            _dir = os.path.join('files', func_name)

            if not os.path.exists(_dir):
                os.makedirs(_dir, exist_ok=True)

            self.path = os.path.join(_dir, f'{file_name}.tif')
            self.json_path = os.path.join(_dir, f'{file_name}.json')

            tifffile.imsave(self.path, data)

            save_data = copy.deepcopy(self.data)
            if self.data.shape[-1] >= 200 and self.data.shape[-2] >= 200:
                save_data = cv2.resize(save_data, (200, 200))

            if len(self.data.shape) == 2:
                self.data = self.data[np.newaxis, :, :]
                save_data = save_data[np.newaxis, :, :]

            images = []
            for i, _img in enumerate(save_data[:10]):
                images.append(_img.tolist())

            pd.DataFrame(images).to_json(self.json_path, indent=4, orient="values")


class TimeSeriesData:
    def __init__(self, data, func_name='timeseries', file_name='timeseries'):
        self.data = data
        _dir = os.path.join('files', func_name)

        if not os.path.exists(_dir):
            os.makedirs(_dir, exist_ok=True)

        if len(self.data.shape) == 1:
            self.data = self.data[np.newaxis, :]

        self.path = os.path.join(_dir, f'{file_name}.json')

        pd.DataFrame(self.data).to_json(self.path, indent=4)


class CorrelationData:
    def __init__(self, data, func_name='heatmap', file_name='heatmap'):
        self.data = data

        _dir = os.path.join('files', func_name)

        if not os.path.exists(_dir):
            os.makedirs(_dir, exist_ok=True)

        self.path = os.path.join(_dir, f'{file_name}.json')

        pd.DataFrame(self.data).to_json(self.path, indent=4, orient="values")


class Suite2pData:
    def __init__(self, data, func_name='suite2p', file_name='suite2p'):
        self.data = data


class IscellData:
    def __init__(self, data, func_name='iscell', file_name='iscell'):
        self.data = data


class RoiData:
    def __init__(self, data, func_name='roi', file_name='roi'):
        self.data = data
