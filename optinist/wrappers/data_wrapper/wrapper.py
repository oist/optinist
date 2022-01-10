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


def get_file_path(func_name, file_name):
    _dir = os.path.join(BASE_DIR, func_name, func_name)
    if not os.path.exists(_dir):
        os.makedirs(_dir, exist_ok=True)

    return os.path.join(_dir, f'{file_name}.json')


def get_images_list(data):
    if len(data.shape) == 2:
        save_data = copy.deepcopy(data)
    elif len(data.shape) == 3:
        save_data = copy.deepcopy(data[:10])
    else:
        assert False, 'data is error'

    if len(data.shape) == 2:
        data = data[np.newaxis, :, :]
        save_data = save_data[np.newaxis, :, :]

    images = []
    for i, _img in enumerate(save_data):
        images.append(_img.tolist())
    
    return images


class ImageData:
    def __init__(self, data, func_name='image', file_name='image'):
        
        if type(data) == str:
            self.path = data
            self.json_path = None
        else:
            _dir = os.path.join(BASE_DIR, func_name)

            # if not os.path.exists(_dir):
            #     os.makedirs(_dir, exist_ok=True)

            # self.json_path = os.path.join(_dir, f'{file_name}.json')
            self.json_path = get_file_path(func_name, file_name)

            self.path = os.path.join(_dir, f'{file_name}.tif')

            tifffile.imsave(self.path, data)

            images = get_images_list(data)

            pd.DataFrame(images).to_json(self.json_path, indent=4, orient="values")

            del images, data
            gc.collect()

    @property
    def data(self):
        return np.array(imageio.volread(self.path))

    def __del__(self):
        del self
        gc.collect()


class TimeSeriesData:
    def __init__(self, data, func_name='timeseries', file_name='timeseries'):
        if type(data) == str:
            self.data = pd.read_csv(data).values.T
        else:
            self.data = data

        if len(self.data.shape) == 1:
            self.data = self.data[np.newaxis, :]

        # if not os.path.exists(_dir):
        #     os.makedirs(_dir, exist_ok=True)

        # self.path = os.path.join(_dir, f'{file_name}.json')
        self.path = get_file_path(func_name, file_name)

        pd.DataFrame(self.data).to_json(self.path, indent=4)

    def __del__(self):
        del self
        gc.collect()


class CorrelationData:
    def __init__(self, data, func_name='heatmap', file_name='heatmap'):
        self.data = data

        # _dir = os.path.join(BASE_DIR, func_name)

        # if not os.path.exists(_dir):
        #     os.makedirs(_dir, exist_ok=True)

        # self.path = os.path.join(_dir, f'{file_name}.json')
        self.path = get_file_path(func_name, file_name)

        pd.DataFrame(self.data).to_json(self.path, indent=4, orient="values")

    def __del__(self):
        del self
        gc.collect()


class RoiData:
    def __init__(self, data, func_name='roi', file_name='roi'):
        self.data = data

        images = get_images_list(data)

        self.path = get_file_path(func_name, file_name)
        pd.DataFrame(images).to_json(self.path, indent=4, orient="values")

        del images, data
        gc.collect()

    def __del__(self):
        del self
        gc.collect()


class Suite2pData:
    def __init__(self, data, func_name='suite2p', file_name='suite2p'):
        self.data = data
        self.path = get_file_path(func_name, file_name)

    def __del__(self):
        del self
        gc.collect()


class IscellData:
    def __init__(self, data, func_name='iscell', file_name='iscell'):
        self.data = data
        self.path = get_file_path(func_name, file_name)

    def __del__(self):
        del self
        gc.collect()
