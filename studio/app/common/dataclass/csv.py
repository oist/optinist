from typing import Optional

import numpy as np
import pandas as pd

from studio.app.common.core.utils.filepath_creater import (
    create_directory,
    join_filepath,
)
from studio.app.common.core.utils.json_writer import JsonWriter
from studio.app.common.dataclass.base import BaseData
from studio.app.common.schemas.outputs import PlotMetaData


class CsvData(BaseData):
    def __init__(
        self, data, params, file_name="csv", meta: Optional[PlotMetaData] = None
    ):
        super().__init__(file_name)
        self.meta = meta

        if isinstance(data, str):
            header = params.get("setHeader", None)
            self.data = pd.read_csv(data, header=header).values

            if params["transpose"]:
                self.data = self.data.T

        else:
            self.data = np.array(data)

        if self.data.ndim == 1:
            self.data = self.data[np.newaxis, :]

    def save_json(self, json_dir):
        # timeseriesだけはdirを返す
        self.json_path = join_filepath([json_dir, self.file_name])
        create_directory(self.json_path)
        JsonWriter.write_plot_meta(json_dir, self.file_name, self.meta)

        for i, data in enumerate(self.data):
            JsonWriter.write_as_split(
                join_filepath([self.json_path, f"{str(i)}.json"]), data
            )
