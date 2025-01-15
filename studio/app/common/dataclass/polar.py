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
            plt.figure()
            ax = plt.subplot(111, polar=True)

            # Convert degrees to radians
            theta_rad = np.radians(self.columns)

            # Basic plot
            ax.plot(theta_rad, self.data[i], color="#636EFA", linewidth=1.5)

            # Set direction and orientation to match plotly defaults
            ax.set_theta_zero_location("E")  # 0° at top
            ax.set_theta_direction(1)  # counterclockwise

            # Set degree labels
            ax.set_xticks(np.radians([0, 45, 90, 135, 180, 225, 270, 315]))
            ax.set_xticklabels(
                ["0°", "45°", "90°", "135°", "180°", "225°", "270°", "315°"]
            )

            # Close the line
            if len(theta_rad) > 0:
                ax.plot(
                    [theta_rad[-1], theta_rad[0]],
                    [self.data[i][-1], self.data[i][0]],
                    color="#636EFA",
                    linewidth=1.5,
                )

            plot_file = join_filepath([output_dir, f"{self.file_name}_{i}.png"])
            plt.tight_layout()
            plt.savefig(plot_file, bbox_inches="tight")
            plt.close()

            save_thumbnail(plot_file)
