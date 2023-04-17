from optinist.api.dataclass.base import BaseData


class HTMLData(BaseData):
    def __init__(self, data, file_name='html'):
        super().__init__(file_name)
        self.data = data

    def set_data_path(self, data_dir):
        return super().set_data_path(data_dir, ext="html")

    def save_data(self):
        with open(self.data_path, "w") as f:
            f.write(self.data)
