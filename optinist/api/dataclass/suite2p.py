from optinist.api.dataclass.base import BaseData
import numpy as np


class Suite2pData(BaseData):
    def __init__(self, data, file_name='suite2p'):
        super().__init__(file_name)
        self.data = data

    def set_data_path(self, data_dir):
        return super().set_data_path(data_dir, ext="npy")

    def save_data(self):
        np.save(self.data_path, self.data)
