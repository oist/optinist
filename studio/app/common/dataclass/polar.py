from typing import Optional

import pandas as pd
import plotly.express as px
import plotly.io as pio

from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.utils.json_writer import JsonWriter
from studio.app.common.core.workflow.workflow import OutputPath, OutputType
from studio.app.common.dataclass.base import BaseData
from studio.app.common.schemas.outputs import PlotMetaData


class PolarData(BaseData):
    def __init__(
        self, data, thetas, file_name="polar", meta: Optional[PlotMetaData] = None
    ):
        # thetas: specify in degrees
        super().__init__(file_name)

        self.data = data
        self.columns = thetas

    def save_json(self, json_dir):
        self.json_path = join_filepath([json_dir, f"{self.file_name}.json"])
        df = pd.DataFrame(self.data, columns=self.columns)
        JsonWriter.write_as_split(self.json_path, df)

    @property
    def output_path(self) -> OutputPath:
        return OutputPath(path=self.json_path, type=OutputType.POLAR)

    def save_plot(self, output_dir):
        for i in range(len(self.data)):
            fig = px.line_polar(
                r=self.data[i],
                theta=self.columns,
                direction="counterclockwise",
                start_angle=0,
                line_close=True,
            )
            pio.write_image(
                fig, join_filepath([output_dir, f"{self.file_name}_{i}.png"])
            )
