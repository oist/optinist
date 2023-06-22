import gc

import imageio
import numpy as np
import tifffile

from optinist.api.dataclass.base import BaseData
from optinist.api.dataclass.utils import create_images_list
from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import create_directory, join_filepath
from optinist.api.utils.json_writer import JsonWriter


class ImageData(BaseData):
    def __init__(self, data, output_dir=DIRPATH.OUTPUT_DIR, file_name="image"):
        super().__init__(file_name)

        self.json_path = None

        if data is None:
            self.path = None
        elif isinstance(data, str):
            self.path = data
        elif isinstance(data, list) and isinstance(data[0], str):
            self.path = data
        else:
            _dir = join_filepath([output_dir, "tiff", file_name])
            create_directory(_dir)

            _path = join_filepath([_dir, f"{file_name}.tif"])
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
        JsonWriter.write_as_split(self.json_path, create_images_list(self.data))


class RoiData(BaseData):
    def __init__(self, data, output_dir=DIRPATH.OUTPUT_DIR, file_name="roi"):
        super().__init__(file_name)

        images = create_images_list(data)

        _dir = join_filepath([output_dir, "tiff", file_name])
        create_directory(_dir)
        self.path = join_filepath([_dir, f"{file_name}.tif"])
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
        JsonWriter.write_as_split(self.json_path, create_images_list(self.data))
