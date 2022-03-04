import os
import numpy as np
import pandas as pd
import imageio
import tifffile
import gc
from cui_api.utils import get_file_path, get_images_list, join_file_path
from cui_api.const import BASE_DIR

class BaseData:
    def __init__(self, file_name):
        self.file_name = file_name

    def save_json(self, json_dir):
        pass


class ImageData(BaseData):
    def __init__(self, data, func_name='image', file_name='image'):
        super().__init__(file_name)

        self.json_path = None

        if type(data) == str:
            self.path = data
        else:
            _dir = join_file_path([BASE_DIR, "tiff", func_name])
            if not os.path.exists(_dir):
                os.makedirs(_dir, exist_ok=True)

            self.path = join_file_path([_dir, f'{file_name}.tif'])

            tifffile.imsave(self.path, data)

            del data
            gc.collect()

    @property
    def data(self):
        return np.array(imageio.volread(self.path))

    def save_json(self, json_dir):
        self.json_path = get_file_path(json_dir, self.file_name)
        images = get_images_list(self.data)
        pd.DataFrame(images).to_json(self.json_path, indent=4, orient="values")
        del images

    def __del__(self):
        del self
        gc.collect()


class TimeSeriesData(BaseData):
    def __init__(self, data, index=None, func_name='timeseries', file_name='timeseries'):
        super().__init__(file_name)

        if type(data) == str:
            self.data = pd.read_csv(data, header=None).values
        else:
            self.data = data

        if len(self.data.shape) == 1:
            self.data = self.data[np.newaxis, :]

        # indexを指定
        self.index = index
        if index == None:
            self.index = np.arange(len(self.data[0]))

    def save_json(self, json_dir):
        # timeseriesだけはdirを返す
        self.json_path = join_file_path([json_dir, self.file_name])
        if not os.path.exists(self.json_path):
            os.makedirs(self.json_path, exist_ok=True)

        for i, data in enumerate(self.data):
            pd.DataFrame(data, index=self.index).to_json(
                join_file_path([self.json_path, f'{str(i)}.json']), indent=4)
            # import pdb; pdb.set_trace()
            # pd.DataFrame(data).to_json(
            #     join_file_path([self.json_path, f'{str(i)}.json']), indent=4)

    def __del__(self):
        del self
        gc.collect()


class CsvData(BaseData):
    def __init__(self, data, params, func_name='csv', file_name='csv'):
        super().__init__(file_name)

        if type(data) == str:
            self.data = pd.read_csv(data, header=None).values

            if params["transpose"]:
                self.data = self.data.T

            if params["setColumn"] is not None:
                header = params["setColumn"]
                self.data = self.data[header:]
        else:
            self.data = data

        if len(self.data.shape) == 1:
            self.data = self.data[np.newaxis, :]

    def save_json(self, json_dir):
        # timeseriesだけはdirを返す
        self.json_path = join_file_path([json_dir, self.file_name])
        if not os.path.exists(self.json_path):
            os.makedirs(self.json_path, exist_ok=True)

        for i, data in enumerate(self.data):
            pd.DataFrame(data).to_json(
                join_file_path([self.json_path, f'{str(i)}.json']), indent=4)

    def __del__(self):
        del self
        gc.collect()


class CorrelationData(BaseData):
    def __init__(self, data, func_name='heatmap', file_name='heatmap'):
        super().__init__(file_name)

        self.data = data

    def save_json(self, json_dir):
        self.json_path = get_file_path(json_dir, self.file_name)
        pd.DataFrame(self.data).to_json(self.json_path, indent=4, orient="values")

    def __del__(self):
        del self
        gc.collect()


class RoiData(BaseData):
    def __init__(self, data, func_name='roi', file_name='roi'):
        super().__init__(file_name)

        images = get_images_list(data)

        _dir = join_file_path([BASE_DIR, "tiff", func_name])
        if not os.path.exists(_dir):
            os.makedirs(_dir, exist_ok=True)
        self.path = join_file_path([_dir, f'{file_name}.tif'])
        tifffile.imsave(self.path, images)

        del images, data
        gc.collect()

    @property
    def data(self):
        data = np.array(imageio.volread(self.path))
        if data.ndim == 3:
            return data[0]
        elif data.ndim == 2:
            return data

    def save_json(self, json_dir):
        self.json_path = get_file_path(json_dir, self.file_name)
        images = get_images_list(self.data)
        pd.DataFrame(images).to_json(self.json_path, indent=4, orient="values")

    def __del__(self):
        del self
        gc.collect()


class Suite2pData(BaseData):
    def __init__(self, data, func_name='suite2p', file_name='suite2p'):
        super().__init__(file_name)

        self.data = data
        _dir = join_file_path([BASE_DIR, func_name])
        self.path = get_file_path(_dir, file_name)

    def __del__(self):
        del self
        gc.collect()


class IscellData(BaseData):
    def __init__(self, data, func_name='iscell', file_name='iscell'):
        super().__init__(file_name)

        self.data = data
        _dir = join_file_path([BASE_DIR, func_name])
        self.path = get_file_path(_dir, file_name)

    def __del__(self):
        del self
        gc.collect()


class ScatterData(BaseData):
    def __init__(self, data, func_name='scatter', file_name='scatter'):
        super().__init__(file_name)

        if not data.ndim == 2:
            raise 'Scatter Dimension Error'

        self.data = data

    def save_json(self, json_dir):
        self.json_path = get_file_path(json_dir, self.file_name)
        pd.DataFrame(self.data).to_json(self.json_path, indent=4)

    def __del__(self):
        del self
        gc.collect()


class BarData(BaseData):
    def __init__(self, data, func_name='bar', file_name='bar'):
        super().__init__(file_name)

        data = np.array(data)

        if not data.ndim == 1:
            raise 'Bar Dimension Error'

        self.data = data

    def save_json(self, json_dir):
        self.json_path = get_file_path(json_dir, self.file_name)
        pd.DataFrame(self.data).to_json(self.json_path, indent=4)

    def __del__(self):
        del self
        gc.collect()
