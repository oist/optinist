from typing import Optional

from studio.app.common.schemas.outputs import PlotMetaData
from studio.app.optinist.dataclass.fluo import FluoData


class SpikingActivityData(FluoData):
    def __init__(
        self,
        data,
        std=None,
        index=None,
        params=None,
        file_name="spiking_activity",
        meta: Optional[PlotMetaData] = None,
    ):
        super().__init__(
            data=data,
            std=std,
            index=index,
            params=params,
            file_name=file_name,
            meta=meta,
        )
