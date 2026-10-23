import json
import mmap
import os
from abc import ABC, abstractmethod
from io import BufferedReader
from pathlib import Path

from typing_extensions import Optional

from studio.app.common.schemas.outputs import (
    JsonTimeSeriesData,
    OutputData,
    PaginatedLineResult,
    PlotMetaData,
    TextPosition,
)


def get_folder_size(folder: str):
    return sum(file.stat().st_size for file in Path(folder).rglob("*"))


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
    def read_as_timeseries_from_df(cls, df) -> JsonTimeSeriesData:
        """
        Create JsonTimeSeriesData from a pandas DataFrame.
        Expects DataFrame with index as xrange and columns
          containing 'data' and optionally 'std'.
        Handles None/null values (which may come from NaN in original data).
        """
        import math

        import pandas as pd

        xrange = [str(idx) for idx in df.index.tolist()]

        # Convert data, handling NaN/None values (null in JSON)
        data = {}
        for idx, val in zip(df.index, df["data"]):
            # Check for NaN/None and convert to None (null in JSON)
            if pd.isna(val):
                data[str(idx)] = None
            elif isinstance(val, float) and (math.isnan(val) or math.isinf(val)):
                data[str(idx)] = None
            else:
                data[str(idx)] = float(val)

        std = None
        if "std" in df.columns:
            std = {}
            for idx, val in zip(df.index, df["std"]):
                # Check for NaN/None and convert to None (null in JSON)
                if pd.isna(val):
                    std[str(idx)] = None
                elif isinstance(val, float) and (math.isnan(val) or math.isinf(val)):
                    std[str(idx)] = None
                else:
                    std[str(idx)] = float(val)

        return JsonTimeSeriesData(
            xrange=xrange,
            data=data,
            std=std,
        )

    @classmethod
    def read_as_plot_meta(cls, filepath) -> PlotMetaData:
        json_data = cls.read(filepath) if os.path.exists(filepath) else {}
        return PlotMetaData(**json_data)


class ContentUnitReader(ABC):
    """Abstract base class for defining how to read content units"""

    @abstractmethod
    def is_unit_start(self, line: bytes) -> bool:
        """Determine if a line starts a new content unit"""
        pass

    @abstractmethod
    def parse(self, content: bytes) -> dict:
        """Parse a content unit into a structured format"""
        pass

    def validate(self, content: bytes) -> bool:
        """
        Condition to check whether line should be included in output data.
        """
        return True


class LineReader(ContentUnitReader):
    """Simple line reader that treats each line as a unit"""

    def is_unit_start(self, line: bytes) -> bool:
        return True

    def parse(self, content: bytes) -> dict:
        return {"raw": content}


class PaginatedFileReader:
    def __init__(self, file_path, **kwargs):
        if not os.path.exists(file_path):
            raise Exception(f"{file_path} does not exist.")
        self.file_path = file_path
        self.unit_reader = LineReader()

    def _read_forward(
        self,
        file: BufferedReader,
        offset: int,
        stop_offset: Optional[int],
        limit: int = 50,
        include_pattern: Optional[bytes] = None,
        search_match_case=True,
    ) -> PaginatedLineResult:
        file.seek(offset)

        units = []
        current_unit_buffer = b""
        while len(units) < (limit if limit != 0 else float("inf")):
            line = file.readline()
            if not line:
                break

            if self.unit_reader.is_unit_start(line):
                if current_unit_buffer:
                    if self.unit_reader.validate(current_unit_buffer):
                        if include_pattern:
                            if (
                                not search_match_case
                                and include_pattern.lower()
                                in current_unit_buffer.lower()
                            ):
                                units.append(current_unit_buffer)
                            elif include_pattern in current_unit_buffer:
                                units.append(current_unit_buffer)
                        else:
                            units.append(current_unit_buffer)
                current_unit_buffer = line
            else:
                current_unit_buffer += line

            if stop_offset and file.tell() > stop_offset:
                break

        next_offset = file.tell() - len(current_unit_buffer)
        data = [line.decode().strip() for line in units]

        return PaginatedLineResult(
            next_offset=next_offset,
            prev_offset=offset,
            data=data,
        )

    def _read_backward(
        self,
        file: BufferedReader,
        offset: int,
        limit: int = 50,
        read_chunk_size=1024,
        include_pattern: Optional[bytes] = None,
        search_match_case=True,
    ) -> PaginatedLineResult:
        file.seek(0, 2)
        file_size = file.tell()
        next_offset = offset = min(offset, file_size)

        segment = b""
        current_unit_buffer = b""
        units = []
        while offset > 0 and len(units) < limit:
            chunk_size = min(read_chunk_size, offset)
            offset -= chunk_size
            file.seek(offset)

            buffer = file.read(chunk_size)

            # Append previous segment
            if segment:
                buffer += segment
                segment = b""

            buffer_lines = buffer.splitlines(keepends=True)

            # first line may be incomplete, save it as segment to append to next chunk
            if buffer_lines and offset > 0:
                segment = buffer_lines.pop(0)

            # Prepend valid lines to `lines`
            len_buf_lines = len(buffer_lines)
            for i in range(len_buf_lines - 1, -1, -1):
                line = buffer_lines.pop(i)

                current_unit_buffer = line + current_unit_buffer
                if self.unit_reader.is_unit_start(line):
                    if self.unit_reader.validate(current_unit_buffer):
                        if include_pattern:
                            if (
                                not search_match_case
                                and include_pattern.lower()
                                in current_unit_buffer.lower()
                            ):
                                units.insert(0, current_unit_buffer)
                            elif include_pattern in current_unit_buffer:
                                units.insert(0, current_unit_buffer)
                        else:
                            units.insert(0, current_unit_buffer)
                        current_unit_buffer = b""
                        if len(units) == limit:
                            segment += b"".join(buffer_lines)
                            break
                    else:
                        current_unit_buffer = b""

        prev_offset = offset + len(segment)

        data = [line.decode().strip() for line in units]

        return PaginatedLineResult(
            next_offset=next_offset,
            prev_offset=prev_offset,
            data=data,
        )

    def read_from_offset(
        self,
        offset: int,
        stop_offset: Optional[int] = None,
        limit: int = 50,
        reverse: bool = False,
    ) -> PaginatedLineResult:
        with open(self.file_path, "rb") as file:
            if offset == -1:
                file.seek(0, 2)
                offset = file.tell()
            if stop_offset:
                limit = 0
                if stop_offset == -1:
                    file.seek(0, 2)
                    stop_offset = file.tell()

            if reverse:
                return self._read_backward(file, offset, limit)
            else:
                return self._read_forward(file, offset, stop_offset, limit)

    def _find_case_sensitive(self, mm: mmap.mmap, search_bytes: bytes, offset, reverse):
        return (
            mm.rfind(search_bytes, 0, offset)
            if reverse
            else mm.find(search_bytes, offset)
        )

    def _find_case_insensitive(
        self, mm: mmap.mmap, search_bytes: bytes, offset, reverse, chunk_size=4096
    ):
        """Find text case-insensitively by processing in chunks."""
        search_len = len(search_bytes)
        file_size = len(mm)
        search_bytes = search_bytes.lower()

        if reverse:
            pos = min(offset, file_size)
        else:
            pos = max(0, offset)

        while (pos >= 0 and reverse) or (pos < file_size and not reverse):
            start_pos = max(0, pos - chunk_size) if reverse else pos
            end_pos = min(pos + chunk_size + search_len, file_size)

            chunk = mm[start_pos:end_pos].lower()
            found_pos = (
                chunk.rfind(search_bytes) if reverse else chunk.find(search_bytes)
            )

            if found_pos != -1:
                return start_pos + found_pos

            if reverse:
                pos -= chunk_size - search_len + 1
                if pos < 0:
                    break
            else:
                pos += chunk_size - search_len + 1

        return -1

    def _find_line_boundaries(self, mm: mmap.mmap, pos, offset, reverse):
        """Find start and end of line for the matched text."""
        start_of_line, end_of_line = None, None

        if reverse:
            eol = mm.find(b"\n", pos, offset)
            end_of_line = offset if eol == -1 else eol
        else:
            sol = mm.rfind(b"\n", 0, pos)
            start_of_line = 0 if sol == -1 else sol

        return start_of_line, end_of_line

    def _get_text_position(
        self,
        search_text: str,
        offset: int,
        reverse: bool = False,
        search_match_case=True,
    ):
        position = TextPosition()
        with open(self.file_path, "r+b") as file:
            if offset == -1:
                file.seek(0, 2)
                offset = file.tell()

            with mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                search_bytes = search_text.encode()

                pos = (
                    self._find_case_sensitive(mm, search_bytes, offset, reverse)
                    if search_match_case
                    else self._find_case_insensitive(mm, search_bytes, offset, reverse)
                )
                if pos == -1:
                    return position, offset
                position.pos = pos
                (
                    position.start_of_line,
                    position.end_of_line,
                ) = self._find_line_boundaries(mm, pos, offset, reverse)

        return (position, offset)

    def get_unit_position_from_search_text(
        self, search: str, offset: int, reverse, search_match_case=True
    ):
        text_position, offset = self._get_text_position(
            search, offset, reverse, search_match_case=False
        )

        if text_position.pos is None:
            return None, offset

        with open(self.file_path, "rb") as file:
            if reverse:
                logs = self._read_backward(
                    file,
                    offset=text_position.end_of_line,
                    limit=1,
                    include_pattern=search.encode(),
                    search_match_case=search_match_case,
                )
                if not logs.data:
                    return None, offset
                return logs.prev_offset, offset
            else:
                logs = self._read_forward(
                    file,
                    offset=text_position.start_of_line,
                    stop_offset=None,
                    limit=1,
                    include_pattern=search.encode(),
                    search_match_case=search_match_case,
                )
                if not logs.data:
                    return None, offset
                return logs.next_offset, offset
