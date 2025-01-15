from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.utils.json_writer import JsonWriter
from studio.app.common.core.workflow.workflow import OutputPath, OutputType
from studio.app.common.dataclass.base import BaseData
from studio.app.common.dataclass.utils import save_thumbnail
from studio.app.common.schemas.outputs import PlotMetaData


class HistogramData(BaseData):
    def __init__(
        self,
        data,
        file_name="histogram",
        meta: Optional[PlotMetaData] = None,
    ):
        super().__init__(file_name)

        if isinstance(data, list):
            data = np.array(data)
        assert isinstance(data, np.ndarray), "Histogram Type Error"
        assert data.ndim == 1, "Histogram Dimension Error"
        self.data = data.reshape(1, -1)

    def save_json(self, json_dir):
        self.json_path = join_filepath([json_dir, f"{self.file_name}.json"])
        df = pd.DataFrame(self.data)
        JsonWriter.write_as_split(self.json_path, df)

    @property
    def output_path(self) -> OutputPath:
        return OutputPath(path=self.json_path, type=OutputType.HISTOGRAM)

    def save_plot(self, output_dir):
        plt.figure()
        plt.hist(self.data[0], bins=20)
        plot_file = join_filepath([output_dir, f"{self.file_name}.png"])
        plt.savefig(plot_file, bbox_inches="tight")
        plt.close()

        save_thumbnail(plot_file)
