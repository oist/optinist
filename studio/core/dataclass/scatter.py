from studio.core.dataclass.base import BaseData
from studio.core.utils.filepath_creater import join_filepath
from studio.core.utils.json_writer import JsonWriter
from studio.core.workflow.workflow import OutputPath, OutputType


class ScatterData(BaseData):
    def __init__(self, data, file_name="scatter"):
        super().__init__(file_name)

        assert data.ndim <= 2, "Scatter Dimension Error"

        self.data = data.T

    def save_json(self, json_dir):
        self.json_path = join_filepath([json_dir, f"{self.file_name}.json"])
        JsonWriter.write_as_split(self.json_path, self.data)

    @property
    def output_path(self) -> OutputPath:
        return OutputPath(path=self.json_path, type=OutputType.SCATTER)
