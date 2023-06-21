from optinist.api.dataclass.base import BaseData
from optinist.api.utils.filepath_creater import join_filepath
from optinist.api.utils.json_writer import JsonWriter


class ScatterData(BaseData):
    def __init__(self, data, file_name="scatter"):
        super().__init__(file_name)

        assert data.ndim <= 2, "Scatter Dimension Error"

        self.data = data.T

    def save_json(self, json_dir):
        self.json_path = join_filepath([json_dir, f"{self.file_name}.json"])
        JsonWriter.write_as_split(self.json_path, self.data)
