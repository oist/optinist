import gc
from typing import Dict, List, Optional

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
from studio.app.common.dataclass.image import ImageData
from studio.app.common.dataclass.utils import create_images_list
from studio.app.common.schemas.outputs import PlotMetaData
from studio.app.dir_path import DIRPATH
from studio.app.optinist.schemas.roi import RoiPos, RoiStatus


class RoiData(BaseData):
    def __init__(
        self,
        data,
        output_dir=DIRPATH.OUTPUT_DIR,
        file_name="roi",
        meta: Optional[PlotMetaData] = None,
    ):
        super().__init__(file_name)
        self.meta = meta

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
        JsonWriter.write_plot_meta(json_dir, self.file_name, self.meta)

    @property
    def output_path(self) -> OutputPath:
        return OutputPath(path=self.json_path, type=OutputType.ROI)


class EditRoiData(BaseData):
    def __init__(self, images, im):
        self.images: ImageData = images
        self.im = im
        self.temp_add_roi: Dict[int, RoiPos] = {}
        self.temp_merge_roi: Dict[float, List[int]] = {}
        self.temp_delete_roi: Dict[float, None] = {}
        self.add_roi = []
        self.merge_roi = []
        self.delete_roi = []

    @property
    def temp_merge_roi_list(self) -> list:
        merge_roi = [(k, *v, -1.0) for k, v in self.temp_merge_roi.items()]
        return [id for ids in merge_roi for id in ids]

    def status(self):
        return RoiStatus(
            temp_add_roi=list(self.temp_add_roi.keys()),
            temp_merge_roi=self.temp_merge_roi_list,
            temp_delete_roi=list(self.temp_delete_roi.keys()),
        )

    def commit(self):
        self.add_roi += self.temp_add_roi.keys()
        self.delete_roi += self.temp_delete_roi.keys()
        self.merge_roi += self.temp_merge_roi_list

        self.temp_add_roi = {}
        self.temp_merge_roi = {}
        self.temp_delete_roi = {}

    def cancel(self):
        self.temp_add_roi = {}
        self.temp_merge_roi = {}
        self.temp_delete_roi = {}
