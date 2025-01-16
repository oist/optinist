from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd

from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.utils.json_writer import JsonWriter
from studio.app.common.core.workflow.workflow import OutputPath, OutputType
from studio.app.common.dataclass.base import BaseData
from studio.app.common.dataclass.utils import save_thumbnail
from studio.app.common.schemas.outputs import PlotMetaData


class LineData(BaseData):
    def __init__(
        self, data, columns, file_name="line", meta: Optional[PlotMetaData] = None
    ):
        super().__init__(file_name)

        self.data = data
        self.columns = columns

    def save_json(self, json_dir):
        self.json_path = join_filepath([json_dir, f"{self.file_name}.json"])
        df = pd.DataFrame(self.data, columns=self.columns)
        JsonWriter.write_as_split(self.json_path, df)

    @property
    def output_path(self) -> OutputPath:
        return OutputPath(path=self.json_path, type=OutputType.LINE)

    def save_plot(self, output_dir):
        for i in range(len(self.data)):
            plt.figure()
            plt.plot(self.columns, self.data[i])
            plt.grid(True, linestyle="--", alpha=0.7)
            plot_file = join_filepath([output_dir, f"{self.file_name}_{i}.png"])
            plt.savefig(plot_file, bbox_inches="tight")
            plt.close()

            save_thumbnail(plot_file)
