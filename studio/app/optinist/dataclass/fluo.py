from typing import Optional

from studio.app.common.dataclass.timeseries import TimeSeriesData
from studio.app.common.schemas.outputs import PlotMetaData


class FluoData(TimeSeriesData):
    def __init__(
        self,
        data,
        std=None,
        index=None,
        params=None,
        file_name="fluo",
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
