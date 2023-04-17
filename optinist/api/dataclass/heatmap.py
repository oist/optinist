import numpy as np
import pandas as pd

from optinist.api.utils.json_writer import JsonWriter

from optinist.api.dataclass.base import BaseData


class HeatMapData(BaseData):
    def __init__(self, data, columns=None, file_name='heatmap'):
        super().__init__(file_name)
        self.data = data

        # indexを指定
        if columns is not None:
            self.columns = columns
        else:
            self.columns = np.arange(len(self.data[0]))

    def save_data(self):
        df = pd.DataFrame(
            self.data,
            columns=self.columns,
        )
        JsonWriter.write_as_split(self.data_path, df)
