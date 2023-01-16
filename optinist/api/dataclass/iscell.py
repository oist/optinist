from optinist.api.dataclass.base import BaseData
import os
import numpy as np

class IscellData(BaseData):
    def __init__(self, data, file_name='iscell'):
        super().__init__(file_name)
        self.data = data

    def save_json(self, json_dir):
        self.json_path = os.path.join(json_dir,f"{self.file_name}.npy")
        np.save(self.json_path, self.data)