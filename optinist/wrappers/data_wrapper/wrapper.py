import os
import numpy as np
import pandas as pd
import imageio
from PIL import Image
import cv2
import tifffile
import copy
import gc

BASE_DIR = '/tmp/optinist'

class ImageData:
    def __init__(self, data, func_name='image', file_name='image'):
        
        if type(data) == str:
            self.path = data
            self.json_path = None
        else:
            _dir = os.path.join(BASE_DIR, func_name)

            if not os.path.exists(_dir):
                os.makedirs(_dir, exist_ok=True)

            self.path = os.path.join(_dir, f'{file_name}.tif')
            self.json_path = os.path.join(_dir, f'{file_name}.json')

            tifffile.imsave(self.path, data)

            if len(data.shape) == 2:
                save_data = copy.deepcopy(data)
            else:
                save_data = copy.deepcopy(data[:10])

            if data.shape[-1] >= 200 and data.shape[-2] >= 200:
                save_data = cv2.resize(save_data, (200, 200))

            if len(data.shape) == 2:
                data = data[np.newaxis, :, :]
                save_data = save_data[np.newaxis, :, :]

            images = []
            for i, _img in enumerate(save_data):
                images.append(_img.tolist())

            pd.DataFrame(images).to_json(self.json_path, indent=4, orient="values")

            del save_data, images, data
            gc.collect()

    @property
    def data(self):
        return np.array(imageio.volread(self.path))

    def __del__(self):
        print('deleting')
        del self
        gc.collect()

class TimeSeriesData:
    def __init__(self, data, func_name='timeseries', file_name='timeseries'):
        self.data = data
        _dir = os.path.join(BASE_DIR, func_name)

        if not os.path.exists(_dir):
            os.makedirs(_dir, exist_ok=True)

        if len(self.data.shape) == 1:
            self.data = self.data[np.newaxis, :]

        self.path = os.path.join(_dir, f'{file_name}.json')

        pd.DataFrame(self.data).to_json(self.path, indent=4)

    def __del__(self):
        print('deleting')
        del self
        gc.collect()

class CorrelationData:
    def __init__(self, data, func_name='heatmap', file_name='heatmap'):
        self.data = data

        _dir = os.path.join(BASE_DIR, func_name)

        if not os.path.exists(_dir):
            os.makedirs(_dir, exist_ok=True)

        self.path = os.path.join(_dir, f'{file_name}.json')

        pd.DataFrame(self.data).to_json(self.path, indent=4, orient="values")

    def __del__(self):
        print('deleting')
        del self
        gc.collect()

class Suite2pData:
    def __init__(self, data, func_name='suite2p', file_name='suite2p'):
        self.data = data

    def __del__(self):
        print('deleting')
        del self
        gc.collect()

class IscellData:
    def __init__(self, data, func_name='iscell', file_name='iscell'):
        self.data = data

    def __del__(self):
        print('deleting')
        del self
        gc.collect()

class RoiData:
    def __init__(self, data, func_name='roi', file_name='roi'):
        self.data = data

    def __del__(self):
        print('deleting')
        del self
        gc.collect()
