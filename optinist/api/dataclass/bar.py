import numpy as np
import pandas as pd

from optinist.api.dataclass.base import BaseData
from optinist.api.utils.filepath_creater import join_filepath
from optinist.api.utils.json_writer import JsonWriter


class BarData(BaseData):
    def __init__(self, data, index=None, file_name="bar"):
        super().__init__(file_name)
        data = np.array(data)

        assert data.ndim <= 2, "Bar Dimension Error"

        if data.ndim == 1:
            data = data[np.newaxis]

        assert data.ndim == 2, "Bar Dimesion is not 2"

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
