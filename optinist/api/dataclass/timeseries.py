import numpy as np
import pandas as pd

from optinist.api.dataclass.base import BaseData
from optinist.api.utils.filepath_creater import create_directory, join_filepath
from optinist.api.utils.json_writer import JsonWriter


class TimeSeriesData(BaseData):
    def __init__(
        self, data, std=None, index=None, cell_numbers=None, file_name="timeseries"
    ):
        super().__init__(file_name)

        assert data.ndim <= 2, "TimeSeries Dimension Error"

        if isinstance(data, str):
            self.data = pd.read_csv(data, header=None).values
        else:
            self.data = data

        if len(self.data.shape) == 1:
            self.data = self.data[np.newaxis, :]

        self.std = std

        # indexを指定
        if index is not None:
            self.index = index
        else:
            self.index = np.arange(len(self.data[0]))

        # cell番号を表示
        if cell_numbers is not None:
            self.cell_numbers = cell_numbers + 1
        else:
            self.cell_numbers = range(1, len(self.data) + 1)

    def save_json(self, json_dir):
        # timeseriesだけはdirを返す
        self.json_path = join_filepath([json_dir, self.file_name])
        create_directory(self.json_path, delete_dir=True)

        for i, cell_i in enumerate(self.cell_numbers):
            data = self.data[i]
            if self.std is not None:
                std = self.std[i]
                df = pd.DataFrame(
                    np.concatenate([data[:, np.newaxis], std[:, np.newaxis]], axis=1),
                    index=self.index,
                    columns=["data", "std"],
                )
            else:
                df = pd.DataFrame(data, index=self.index, columns=["data"])

            JsonWriter.write(join_filepath([self.json_path, f"{str(cell_i)}.json"]), df)


class FluoData(TimeSeriesData):
    def __init__(self, data, std=None, index=None, file_name="fluo"):
        super().__init__(
            data=data,
            std=std,
            index=index,
            file_name=file_name,
        )


class BehaviorData(TimeSeriesData):
    def __init__(self, data, std=None, index=None, file_name="behavior"):
        super().__init__(
            data=data,
            std=std,
            index=index,
            file_name=file_name,
        )
