import numpy as np
import pandas as pd

from studio.api.utils.filepath_creater import join_filepath
from studio.api.utils.json_writer import JsonWriter

from studio.api.dataclass.base import BaseData


class HeatMapData(BaseData):
    def __init__(self, data, columns=None, file_name='heatmap'):
        super().__init__(file_name)
        self.data = data

        # indexを指定
        if columns is not None:
            self.columns = columns
        else:
            self.columns = np.arange(len(self.data[0]))

    def save_json(self, json_dir):
        self.json_path = join_filepath([json_dir, f"{self.file_name}.json"])
        df = pd.DataFrame(
            self.data,
            columns=self.columns,
        )
        JsonWriter.write_as_split(self.json_path, df)
