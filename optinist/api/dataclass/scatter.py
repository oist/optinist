from optinist.api.dataclass.base import BaseData
from optinist.api.utils.json_writer import JsonWriter


class ScatterData(BaseData):
    def __init__(self, data, file_name='scatter'):
        super().__init__(file_name)

        assert data.ndim <= 2, 'Scatter Dimension Error'

        self.data = data.T

    def save_data(self):
        JsonWriter.write_as_split(self.data_path, self.data)
