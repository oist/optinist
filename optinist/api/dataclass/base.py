import gc
from optinist.api.utils.filepath_creater import join_filepath


class BaseData:
    def __init__(self, file_name):
        self.file_name = file_name
        self.data_path = None

    def set_data_path(self, data_dir, ext="json"):
        self.data_path = join_filepath([data_dir, f"{self.file_name}.{ext}"])

    def save_data(self):
        pass

    def __del__(self):
        del self
        gc.collect()
