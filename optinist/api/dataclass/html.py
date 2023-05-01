from optinist.api.dataclass.base import BaseData
from optinist.api.utils.filepath_creater import join_filepath


class HTMLData(BaseData):
    def __init__(self, data, file_name="html"):
        super().__init__(file_name)
        self.data = data

    def save_json(self, json_dir):
        self.json_path = join_filepath([json_dir, f"{self.file_name}.html"])

        with open(self.json_path, "w") as f:
            f.write(self.data)
