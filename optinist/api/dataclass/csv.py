import numpy as np
import pandas as pd

from optinist.api.dataclass.base import BaseData
from optinist.api.utils.filepath_creater import create_directory, join_filepath
from optinist.api.utils.json_writer import JsonWriter


class CsvData(BaseData):
    def __init__(self, data, params, file_name="csv"):
        super().__init__(file_name)

        if isinstance(data, str):
            self.data = pd.read_csv(data, header=None).values

            if params["transpose"]:
                self.data = self.data.T

            if params["setHeader"] is not None:
                header = params["setHeader"]
                self.data = self.data[header:]
        else:
            self.data = np.array(data)

        if self.data.ndim == 1:
            self.data = self.data[np.newaxis, :]

    def save_json(self, json_dir):
        # timeseriesだけはdirを返す
        self.json_path = join_filepath([json_dir, self.file_name])
        create_directory(self.json_path)

        for i, data in enumerate(self.data):
            JsonWriter.write_as_split(
                join_filepath([self.json_path, f"{str(i)}.json"]), data
            )
