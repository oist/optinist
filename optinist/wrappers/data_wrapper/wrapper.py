import os
import numpy as np
import pandas as pd
import imageio
from PIL import Image
import cv2
import tifffile
import gc
from .utils import get_file_path, get_images_list

BASE_DIR = '/tmp/optinist'


class ImageData:
    def __init__(self, data, func_name='image', file_name='image'):
        
        if type(data) == str:
            self.path = data
            self.json_path = None
        else:
            _dir = os.path.join(BASE_DIR, func_name)

            self.json_path = get_file_path(_dir, file_name)

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
            self.data = pd.read_csv(data).values
        else:
            self.data = data

        if len(self.data.shape) == 1:
            self.data = self.data[np.newaxis, :]

        _dir = os.path.join(BASE_DIR, func_name, f'{file_name}')
        # self.path = get_file_path(_dir, file_name)
        if not os.path.exists(_dir):
            os.makedirs(_dir, exist_ok=True)

        self.path = os.path.join(_dir)
        # import pdb; pdb.set_trace()
        pd.DataFrame(self.data.T).to_csv(os.path.join(self.path, f'{file_name}.csv'))

        for i, data in enumerate(self.data):
            pd.DataFrame(data).to_json(
                os.path.join(_dir, f'{str(i)}.json'), indent=4)

    def __del__(self):
        del self
        gc.collect()


class CorrelationData:
    def __init__(self, data, func_name='heatmap', file_name='heatmap'):
        self.data = data

        _dir = os.path.join(BASE_DIR, func_name)
        self.path = get_file_path(_dir, file_name)

        pd.DataFrame(self.data).to_json(self.path, indent=4, orient="values")

    def __del__(self):
        del self
        gc.collect()


class RoiData:
    def __init__(self, data, func_name='roi', file_name='roi'):
        self.data = data

        images = get_images_list(data)

        _dir = os.path.join(BASE_DIR, func_name)
        self.path = get_file_path(_dir, file_name)
        pd.DataFrame(images).to_json(self.path, indent=4, orient="values")

        del images, data
        gc.collect()

    def __del__(self):
        del self
        gc.collect()


class Suite2pData:
    def __init__(self, data, func_name='suite2p', file_name='suite2p'):
        self.data = data
        _dir = os.path.join(BASE_DIR, func_name)
        self.path = get_file_path(_dir, file_name)

    def __del__(self):
        del self
        gc.collect()


class IscellData:
    def __init__(self, data, func_name='iscell', file_name='iscell'):
        self.data = data
        _dir = os.path.join(BASE_DIR, func_name)
        self.path = get_file_path(_dir, file_name)

    def __del__(self):
        del self
        gc.collect()


class ScatterData:
    def __init__(self, data, func_name='scatter', file_name='scatter'):

        if not (data.ndim == 2 and data.shape[1] == 2):
            raise 'Scatter Dimension Error'

        _dir = os.path.join(BASE_DIR, func_name, f'{file_name}')
        if not os.path.exists(_dir):
            os.makedirs(_dir, exist_ok=True)

        self.path = get_file_path(_dir, file_name)

        pd.DataFrame(data).to_json(
            os.path.join(_dir, f'{file_name}.json'), indent=4)

    def __del__(self):
        del self
        gc.collect()

