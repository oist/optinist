from optinist.api.dataclass.base import BaseData


class Suite2pData(BaseData):
    def __init__(self, data, file_name='suite2p'):
        super().__init__(file_name)
        self.data = data
