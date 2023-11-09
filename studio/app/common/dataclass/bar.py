from typing import Optional

import numpy as np
import pandas as pd

from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.utils.json_writer import JsonWriter
from studio.app.common.core.workflow.workflow import OutputPath, OutputType
from studio.app.common.dataclass.base import BaseData
from studio.app.common.schemas.outputs import PlotMetaData


class BarData(BaseData):
    def __init__(
        self, data, index=None, file_name="bar", meta: Optional[PlotMetaData] = None
    ):
        super().__init__(file_name)
        self.meta = meta
        data = np.array(data)

        assert data.ndim <= 2, "Bar Dimension Error"

        if data.ndim == 1:
            data = data[np.newaxis]

        assert data.ndim == 2, "Bar Dimension is not 2"

        self.data = data

        # indexを指定
        if index is not None:
            self.index = index
        else:
            self.index = np.arange(len(self.data))

    def save_json(self, json_dir):
        self.json_path = join_filepath([json_dir, f"{self.file_name}.json"])
        df = pd.DataFrame(
            self.data,
            index=self.index,
        )
        JsonWriter.write_as_split(self.json_path, df)
        JsonWriter.write_plot_meta(json_dir, self.file_name, self.meta)

    @property
    def output_path(self) -> OutputPath:
        return OutputPath(path=self.json_path, type=OutputType.BAR)
