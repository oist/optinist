import os
from typing import Optional

import numpy as np

from studio.app.common.core.utils.json_writer import JsonWriter
from studio.app.common.dataclass.base import BaseData
from studio.app.common.schemas.outputs import PlotMetaData


class LccdData(BaseData):
    def __init__(self, data, file_name="lccd", meta: Optional[PlotMetaData] = None):
        super().__init__(file_name)
        self.data = data
        self.meta = meta

    def save_json(self, json_dir):
        self.json_path = os.path.join(json_dir, f"{self.file_name}.npy")
        np.save(self.json_path, self.data)
        JsonWriter.write_plot_meta(json_dir, self.file_name, self.meta)
