import os
import shutil
import numpy as np
import pandas as pd
import imageio
import tifffile
import gc
from pynwb import NWBFile


from optinist.api.utils.filepath_creater import create_directory, join_filepath
from optinist.api.dir_path import DIRPATH
from optinist.api.dataclass.utils import create_images_list
from optinist.api.utils.json_writer import JsonWriter


class NWBFile(NWBFile):
    pass


class BaseData:
    def __init__(self, file_name):
        self.file_name = file_name

    def save_json(self, json_dir):
        pass


class ImageData(BaseData):
    def __init__(self, data, file_name='image'):
        super().__init__(file_name)

        self.json_path = None

        if data is None:
            self.path = None
        elif isinstance(data, str):
            self.path = data
        elif isinstance(data, list) and isinstance(data[0], str):
            self.path = data
        else:
            _dir = join_filepath([DIRPATH.OUTPUT_DIR, "tiff", file_name])
            create_directory(_dir)

            _path = join_filepath([_dir, f'{file_name}.tif'])
            tifffile.imsave(_path, data)
            self.path = [_path]

            del data
            gc.collect()

    @property
    def data(self):
        if isinstance(self.path, list):
            return np.concatenate([imageio.volread(p) for p in self.path])
        else:
            return np.array(imageio.volread(self.path))

    def save_json(self, json_dir):
        self.json_path = join_filepath([json_dir, f"{self.file_name}.json"])
        JsonWriter.write_as_values(
            self.json_path,
            create_images_list(self.data)
        )

    def __del__(self):
        del self
        gc.collect()


class TimeSeriesData(BaseData):
    def __init__(self, data, std=None, index=None, cell_numbers=None, file_name='timeseries'):
        super().__init__(file_name)

        if isinstance(data, str):
            self.data = pd.read_csv(data, header=None).values
        else:
            self.data = data

        if len(self.data.shape) == 1:
            self.data = self.data[np.newaxis, :]

        self.std = std

        # indexを指定        
        if index is not None:
            self.index = index
        else:
            self.index = np.arange(len(self.data[0]))

        # cell番号を表示
        if cell_numbers is not None:
            self.cell_numbers = cell_numbers + 1
        else:
            self.cell_numbers = range(1, len(self.data) + 1)

    def save_json(self, json_dir):
        # timeseriesだけはdirを返す
        self.json_path = join_filepath([json_dir, self.file_name])
        if os.path.exists(self.json_path):
            shutil.rmtree(self.json_path)
        create_directory(self.json_path)

        for i, cell_i in enumerate(self.cell_numbers):
            data = self.data[i]
            if self.std is not None:
                std = self.std[i]
                df = pd.DataFrame(
                    np.concatenate([data[:, np.newaxis], std[:, np.newaxis]], axis=1),
                    index=self.index,
                    columns=["data", "std"]
                )
            else:
                df = pd.DataFrame(
                    data,
                    index=self.index,
                    columns=["data"]
                )

            JsonWriter.write(
                join_filepath([self.json_path, f'{str(cell_i)}.json']),
                df
            )

    def __del__(self):
        del self
        gc.collect()


class FluoData(TimeSeriesData):
    def __init__(self, data, std=None, index=None, file_name='fluo'):
        super().__init__(
            data=data,
            std=std,
            index=index,
            file_name=file_name,
        )


class BehaviorData(TimeSeriesData):
    def __init__(self, data, std=None, index=None, file_name='behavior'):
        super().__init__(
            data=data,
            std=std,
            index=index,
            file_name=file_name,
        )


class CsvData(BaseData):
    def __init__(self, data, params, file_name='csv'):
        super().__init__(file_name)

        if isinstance(data, str):
            self.data = pd.read_csv(data, header=None).values

            if params["transpose"]:
                self.data = self.data.T

            if params["setHeader"] is not None:
                header = params["setHeader"]
                self.data = self.data[header:]
        else:
            self.data = data

        if len(self.data.shape) == 1:
            self.data = self.data[np.newaxis, :]

    def save_json(self, json_dir):
        # timeseriesだけはdirを返す
        self.json_path = join_filepath([json_dir, self.file_name])
        if not os.path.exists(self.json_path):
            os.makedirs(self.json_path)

        for i, data in enumerate(self.data):
            JsonWriter.write(
                join_filepath([self.json_path, f'{str(i)}.json']),
                data
            )

    def __del__(self):
        del self
        gc.collect()


class HeatMapData(BaseData):
    def __init__(self, data, file_name='heatmap'):
        super().__init__(file_name)
        self.data = data

    def save_json(self, json_dir):
        self.json_path = join_filepath([json_dir, f"{self.file_name}.json"])
        JsonWriter.write_as_values(self.json_path, self.data)

    def __del__(self):
        del self
        gc.collect()


class RoiData(BaseData):
    def __init__(self, data, file_name='roi'):
        super().__init__(file_name)

        images = create_images_list(data)

        _dir = join_filepath([DIRPATH.OUTPUT_DIR, "tiff", file_name])
        create_directory(_dir)
        self.path = join_filepath([_dir, f'{file_name}.tif'])
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
        self.json_path = join_filepath([json_dir, f"{self.file_name}.json"])
        JsonWriter.write_as_values(self.json_path, create_images_list(self.data))

    def __del__(self):
        del self
        gc.collect()


class Suite2pData(BaseData):
    def __init__(self, data, file_name='suite2p'):
        super().__init__(file_name)
        self.data = data

    def __del__(self):
        del self
        gc.collect()


class IscellData(BaseData):
    def __init__(self, data, file_name='iscell'):
        super().__init__(file_name)
        self.data = data

    def __del__(self):
        del self
        gc.collect()


class ScatterData(BaseData):
    def __init__(self, data, file_name='scatter'):
        super().__init__(file_name)
        if not data.ndim == 2:
            raise 'Scatter Dimension Error'
        self.data = data

    def save_json(self, json_dir):
        self.json_path = join_filepath([json_dir, f"{self.file_name}.json"])
        JsonWriter.write(self.json_path, self.data)

    def __del__(self):
        del self
        gc.collect()


class BarData(BaseData):
    def __init__(self, data, file_name='bar'):
        super().__init__(file_name)
        data = np.array(data)
        if not data.ndim == 1:
            raise 'Bar Dimension Error'

        self.data = data

    def save_json(self, json_dir):
        self.json_path = join_filepath([json_dir, f"{self.file_name}.json"])
        JsonWriter.write(self.json_path, self.data)

    def __del__(self):
        del self
        gc.collect()


class HTMLData(BaseData):
    def __init__(self, data, file_name='html'):
        super().__init__(file_name)
        self.data = data

    def save_json(self, json_dir):
        self.json_path = join_filepath(json_dir, f"{self.file_name}.html")

        with open(self.json_path, "w") as f:
            f.write(self.data)

    def __del__(self):
        del self
        gc.collect()
