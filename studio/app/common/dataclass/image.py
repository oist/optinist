import gc

import imageio
import numpy as np
import tifffile

from studio.app.common.core.utils.filepath_creater import (
    create_directory,
    join_filepath,
)
from studio.app.common.core.utils.json_writer import JsonWriter
from studio.app.common.core.workflow.workflow import OutputPath, OutputType
from studio.app.common.dataclass.base import BaseData
from studio.app.common.dataclass.utils import create_images_list
from studio.config.dir_path import DIRPATH


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

    @property
    def output_path(self) -> OutputPath:
        return OutputPath(
            path=self.json_path,
            type=OutputType.IMAGE,
            max_index=len(self.data) if self.data.ndim == 3 else 1,
        )
