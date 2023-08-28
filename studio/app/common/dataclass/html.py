from typing import Optional

from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.utils.json_writer import JsonWriter
from studio.app.common.core.workflow.workflow import OutputPath, OutputType
from studio.app.common.dataclass.base import BaseData
from studio.app.common.schemas.outputs import PlotMetaData


class HTMLData(BaseData):
    def __init__(self, data, file_name="html", meta: Optional[PlotMetaData] = None):
        super().__init__(file_name)
        self.data = data
        self.meta = meta

    def save_json(self, json_dir):
        self.json_path = join_filepath([json_dir, f"{self.file_name}.html"])
        JsonWriter.write_plot_meta(json_dir, self.file_name, self.meta)

        with open(self.json_path, "w") as f:
            f.write(self.data)

    @property
    def output_path(self) -> OutputPath:
        return OutputPath(path=self.json_path, type=OutputType.HTML)
