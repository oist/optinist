import gc


class BaseData:
    def __init__(self, file_name):
        self.file_name = file_name

    def save_json(self, json_dir):
        pass

    def __del__(self):
        del self
        gc.collect()
