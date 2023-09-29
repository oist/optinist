from typing import Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio

from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.utils.json_writer import JsonWriter
from studio.app.common.core.workflow.workflow import OutputPath, OutputType
from studio.app.common.dataclass.base import BaseData
from studio.app.common.schemas.outputs import PlotMetaData


class PieData(BaseData):
    def __init__(
        self, data, labels: list, file_name="pie", meta: Optional[PlotMetaData] = None
    ):
        super().__init__(file_name)

        if isinstance(data, list):
            data = np.array(data)
        assert isinstance(data, np.ndarray), "Pie Type Error"
        assert data.ndim == 1, "Pie Dimension Error"
        assert data.shape[0] == len(
            labels
        ), f"labels length is not same as data shape {data.shape}"
        self.data = data.reshape(1, -1)
        self.columns = labels

    def save_json(self, json_dir):
        self.json_path = join_filepath([json_dir, f"{self.file_name}.json"])
        df = pd.DataFrame(self.data, columns=self.columns)
        JsonWriter.write_as_split(self.json_path, df)

    @property
    def output_path(self) -> OutputPath:
        return OutputPath(path=self.json_path, type=OutputType.PIE)

    def save_plot(self, output_dir):
        fig = go.Figure(
            data=[
                go.Pie(
                    labels=self.columns,
                    values=self.data[0],
                    sort=False,
                    direction="clockwise",
                )
            ]
        )
        pio.write_image(fig, join_filepath([output_dir, f"{self.file_name}.png"]))
