import gc

import imageio
import numpy as np
import tifffile

from studio.config.dir_path import DIRPATH
from studio.core.dataclass.base import BaseData
from studio.core.dataclass.utils import create_images_list
from studio.core.utils.filepath_creater import create_directory, join_filepath
from studio.core.utils.json_writer import JsonWriter
from studio.core.workflow.workflow import OutputPath, OutputType


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

    @property
    def output_path(self) -> OutputPath:
        return OutputPath(path=self.json_path, type=OutputType.ROI)
