from studio.app.common.dataclass.base import BaseData


class IscellData(BaseData):
    def __init__(self, data, file_name="iscell"):
        super().__init__(file_name)
        self.data = data
