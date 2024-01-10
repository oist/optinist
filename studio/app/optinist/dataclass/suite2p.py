from typing import Optional

from studio.app.common.dataclass.base import BaseData
from studio.app.common.schemas.outputs import PlotMetaData


class Suite2pData(BaseData):
    def __init__(self, data, file_name="suite2p", meta: Optional[PlotMetaData] = None):
        super().__init__(file_name)
        self.data = data
        self.meta = meta

    def save_json(self, json_dir):
        pass
