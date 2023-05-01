from optinist.api.dataclass.base import BaseData
import os
import numpy as np

class IscellData(BaseData):
    def __init__(self, data, file_name='iscell'):
        super().__init__(file_name)
        self.data = data
