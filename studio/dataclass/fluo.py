from studio.dataclass.timeseries import TimeSeriesData


class FluoData(TimeSeriesData):
    def __init__(self, data, std=None, index=None, file_name="fluo"):
        super().__init__(
            data=data,
            std=std,
            index=index,
            file_name=file_name,
        )
