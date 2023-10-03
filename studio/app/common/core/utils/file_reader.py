import json
import os

from studio.app.common.schemas.outputs import (
    JsonTimeSeriesData,
    OutputData,
    PlotMetaData,
)


class Reader:
    @classmethod
    def read(cls, filepath):
        with open(filepath, "r") as f:
            data = f.read()
        return data

    @classmethod
    def read_as_output(cls, filepath) -> OutputData:
        return OutputData(cls.read(filepath))


class JsonReader:
    @classmethod
    def read(cls, filepath):
        with open(filepath, "r") as f:
            json_data = json.load(f)
        return json_data

    @classmethod
    def read_as_output(cls, filepath) -> OutputData:
        json_data = cls.read(filepath)
        plot_metadata_path = f"{os.path.splitext(filepath)[0]}.plot-meta.json"
        plot_metadata = cls.read_as_plot_meta(plot_metadata_path)

        return OutputData(
            data=json_data["data"],
            columns=json_data["columns"],
            index=json_data["index"],
            meta=plot_metadata,
        )

    @classmethod
    def read_as_timeseries(cls, filepath) -> JsonTimeSeriesData:
        json_data = cls.read(filepath)
        return JsonTimeSeriesData(
            xrange=list(json_data["data"].keys()),
            data=json_data["data"],
            std=json_data["std"] if "std" in json_data else None,
        )

    @classmethod
    def read_as_plot_meta(cls, filepath) -> PlotMetaData:
        json_data = cls.read(filepath) if os.path.exists(filepath) else {}
        return PlotMetaData(**json_data)
