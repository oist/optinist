import os
from typing import Optional

import pandas as pd
from fastapi import APIRouter

from studio.app.common.core.utils.file_reader import JsonReader, Reader
from studio.app.common.core.utils.filepath_creater import (
    create_directory,
    join_filepath,
)
from studio.app.common.core.utils.json_writer import JsonWriter, save_tiff2json
from studio.app.common.dataclass.timeseries_chunk_handler import TimeSeriesChunkHandler
from studio.app.common.schemas.outputs import JsonTimeSeriesData, OutputData
from studio.app.const import ACCEPT_FILE_EXT, ORIGINAL_DATA_EXT
from studio.app.dir_path import DIRPATH

router = APIRouter(prefix="/outputs", tags=["outputs"])


def get_initial_timeseries_data(dirpath) -> JsonTimeSeriesData:
    plot_meta_path = f"{dirpath}.plot-meta.json"
    plot_meta = JsonReader.read_as_plot_meta(plot_meta_path)

    return JsonTimeSeriesData(
        xrange=[],
        data={},
        std={},
        meta=plot_meta,
    )


def _load_timeseries_record(dirpath: str, record_id: str) -> JsonTimeSeriesData:
    """
    Load a single timeseries record from either chunked or legacy format.

    Args:
        dirpath: Directory containing the timeseries data
        record_id: Record identifier (as string)

    Returns:
        JsonTimeSeriesData for the specified record
    """
    if TimeSeriesChunkHandler.is_chunked_format(dirpath):
        # Chunked format
        cell_data = TimeSeriesChunkHandler.get_record_data(dirpath, record_id)
        # Convert from split format to DataFrame
        df = pd.DataFrame(
            cell_data["data"], index=cell_data["index"], columns=cell_data["columns"]
        )
        return JsonReader.read_as_timeseries_from_df(df)
    else:
        # Legacy format
        return JsonReader.read_as_timeseries(
            join_filepath([dirpath, f"{record_id}.json"])
        )


@router.get("/inittimedata/{dirpath:path}", response_model=JsonTimeSeriesData)
async def get_inittimedata(
    dirpath: str,
    isFull: Optional[bool] = None,
):
    full_json_dirpath = dirpath + ORIGINAL_DATA_EXT
    if isFull and os.path.exists(full_json_dirpath):
        dirpath = full_json_dirpath

    # Get all cell indices (supports both chunked and legacy formats)
    file_numbers = TimeSeriesChunkHandler.get_all_record_ids(dirpath)

    # Handle empty case
    if not file_numbers:
        return_data = get_initial_timeseries_data(dirpath)
        return_data.meta = {"title": "0 ROIs found"}  # Set informative message
        return return_data

    # Get first cell data
    index = file_numbers[0]
    str_index = str(index)

    # Load first record using common helper
    json_data = _load_timeseries_record(dirpath, str_index)

    data = {
        str(i): {json_data.xrange[0]: json_data.data[json_data.xrange[0]]}
        for i in file_numbers
    }

    if json_data.std is not None:
        std = {
            str(i): {json_data.xrange[0]: json_data.data[json_data.xrange[0]]}
            for i in file_numbers
        }

    return_data = get_initial_timeseries_data(dirpath)
    return_data.xrange = json_data.xrange
    if json_data.std is not None:
        return_data.std = std

    return_data.data = data
    return_data.data[str_index] = json_data.data
    if json_data.std is not None:
        return_data.std[str_index] = json_data.std

    return return_data


@router.get("/timedata/{dirpath:path}", response_model=JsonTimeSeriesData)
async def get_timedata(
    dirpath: str,
    index: int,
    isFull: Optional[bool] = None,
):
    full_json_dirpath = dirpath + ORIGINAL_DATA_EXT
    if isFull and os.path.exists(full_json_dirpath):
        dirpath = full_json_dirpath

    str_index = str(index)

    # Load record using common helper
    json_data = _load_timeseries_record(dirpath, str_index)

    return_data = get_initial_timeseries_data(dirpath)

    return_data.data[str_index] = json_data.data
    if json_data.std is not None:
        return_data.std[str_index] = json_data.std

    return return_data


@router.get("/alltimedata/{dirpath:path}", response_model=JsonTimeSeriesData)
async def get_alltimedata(dirpath: str):
    return_data = get_initial_timeseries_data(dirpath)

    if TimeSeriesChunkHandler.is_chunked_format(dirpath):
        # Chunked format: load all chunks
        all_records = TimeSeriesChunkHandler.load_all_records(dirpath)

        for cell_index, cell_data in all_records.items():
            # Convert from split format to timeseries format
            df = pd.DataFrame(
                cell_data["data"],
                index=cell_data["index"],
                columns=cell_data["columns"],
            )
            json_data = JsonReader.read_as_timeseries_from_df(df)

            if not return_data.xrange:
                return_data.xrange = json_data.xrange

            return_data.data[cell_index] = json_data.data
            if json_data.std is not None:
                if not return_data.std:
                    return_data.std = {}
                return_data.std[cell_index] = json_data.std
    else:
        # Legacy format: individual files
        from glob import glob

        metadata_files = [
            TimeSeriesChunkHandler.INDEX_MAP_FILENAME,  # chunk_index_map.json
            f"{os.path.basename(dirpath)}.plot-meta.json",
        ]
        for i, path in enumerate(glob(join_filepath([dirpath, "*.json"]))):
            filename = os.path.basename(path)
            # Skip metadata files
            if filename in metadata_files:
                continue

            str_idx = str(os.path.splitext(filename)[0])
            json_data = JsonReader.read_as_timeseries(path)
            if i == 0:
                return_data.xrange = json_data.xrange

            return_data.data[str_idx] = json_data.data
            if json_data.std is not None:
                return_data.std[str_idx] = json_data.std

    return return_data


@router.get("/data/{filepath:path}", response_model=OutputData)
async def get_file(filepath: str):
    return JsonReader.read_as_output(filepath)


@router.get("/html/{filepath:path}", response_model=OutputData)
async def get_html(filepath: str):
    return Reader.read_as_output(filepath)


@router.get("/image/{filepath:path}", response_model=OutputData)
async def get_image(
    filepath: str,
    workspace_id: str,
    start_index: Optional[int] = 0,
    end_index: Optional[int] = 10,
    isFull: Optional[bool] = None,
):
    filename, ext = os.path.splitext(os.path.basename(filepath))

    if filename == "cell_roi" and isFull:
        full_cell_roi_filepath = filepath + ORIGINAL_DATA_EXT
        if os.path.exists(full_cell_roi_filepath):
            filepath = full_cell_roi_filepath

    if ext in ACCEPT_FILE_EXT.TIFF_EXT.value:
        if not filepath.startswith(join_filepath([DIRPATH.OUTPUT_DIR, workspace_id])):
            filepath = join_filepath([DIRPATH.INPUT_DIR, workspace_id, filepath])

        save_dirpath = join_filepath(
            [
                os.path.dirname(filepath),
                filename,
            ]
        )
        json_filepath = join_filepath(
            [save_dirpath, f"{filename}_{str(start_index)}_{str(end_index)}.json"]
        )
        if not os.path.exists(json_filepath):
            save_tiff2json(filepath, save_dirpath, start_index, end_index)
    else:
        json_filepath = filepath

    return JsonReader.read_as_output(json_filepath)


@router.get("/csv/{filepath:path}", response_model=OutputData)
async def get_csv(filepath: str, workspace_id: str):
    filepath = join_filepath([DIRPATH.INPUT_DIR, workspace_id, filepath])

    filename, _ = os.path.splitext(os.path.basename(filepath))
    save_dirpath = join_filepath([os.path.dirname(filepath), filename])
    create_directory(save_dirpath)
    json_filepath = join_filepath([save_dirpath, f"{filename}.json"])

    JsonWriter.write_as_split(json_filepath, pd.read_csv(filepath, header=None))
    return JsonReader.read_as_output(json_filepath)
