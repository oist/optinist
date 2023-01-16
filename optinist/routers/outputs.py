import os
import pandas as pd
from glob import glob
from typing import Optional
from fastapi import APIRouter, status, HTTPException
from optinist.api.dir_path import DIRPATH

from optinist.api.utils.json_writer import JsonWriter, save_tiff2json
from optinist.api.utils.filepath_creater import create_directory, join_filepath
from optinist.routers.const import ACCEPT_TIFF_EXT
from optinist.routers.fileIO.file_reader import JsonReader, Reader
from optinist.routers.model import JsonTimeSeriesData, RoiList, RoiPos, EditRoiSuccess
from optinist.wrappers.suite2p_wrapper.edit_roi import execute_add_ROI, execute_merge_roi, excute_delete_roi

router = APIRouter()


@router.get("/outputs/inittimedata/{dirpath:path}")
async def get_inittimedata(dirpath: str):
    file_numbers = sorted([
        os.path.splitext(os.path.basename(x))[0]
        for x in glob(join_filepath([dirpath, '*.json']))
    ])

    index = file_numbers[0]
    str_index = str(index)

    json_data = JsonReader.read_as_timeseries(
        join_filepath([dirpath, f'{str(index)}.json'])
    )

    return_data = JsonTimeSeriesData(
        xrange=[],
        data={},
        std={},
    )

    data = {
        str(i): {
            json_data.xrange[0]: json_data.data[json_data.xrange[0]]
        }
        for i in file_numbers
    }

    if json_data.std is not None:
        std = {
            str(i): {
                json_data.xrange[0]: json_data.data[json_data.xrange[0]]
            }
            for i in file_numbers
        }

    return_data = JsonTimeSeriesData(
        xrange=json_data.xrange,
        data=data,
        std=std if json_data.std is not None else {},
    )

    return_data.data[str_index] = json_data.data
    if json_data.std is not None:
        return_data.std[str_index] = json_data.std

    return return_data


@router.get("/outputs/timedata/{dirpath:path}")
async def get_timedata(dirpath: str, index: int):
    json_data = JsonReader.read_as_timeseries(
        join_filepath([
            dirpath,
            f'{str(index)}.json'
        ])
    )

    return_data = JsonTimeSeriesData(
        xrange=[],
        data={},
        std={},
    )

    str_index = str(index)
    return_data.data[str_index] = json_data.data
    if json_data.std is not None:
        return_data.std[str_index] = json_data.std

    return return_data


@router.get("/outputs/alltimedata/{dirpath:path}")
async def get_alltimedata(dirpath: str):
    return_data = JsonTimeSeriesData(
        xrange=[],
        data={},
        std={},
    )

    for i, path in enumerate(glob(join_filepath([dirpath, '*.json']))):
        str_idx = str(os.path.splitext(os.path.basename(path))[0])
        json_data = JsonReader.read_as_timeseries(path)
        if i == 0:
            return_data.xrange = json_data.xrange

        return_data.data[str_idx] = json_data.data
        if json_data.std is not None:
            return_data.std[str_idx] = json_data.std

    return return_data


@router.get("/outputs/data/{filepath:path}")
async def get_file(filepath: str):
    return JsonReader.read_as_output(filepath)


@router.get("/outputs/html/{filepath:path}")
async def get_html(filepath: str):
    return Reader.read_as_output(filepath)


@router.get("/outputs/image/{filepath:path}")
async def get_image(
        filepath: str,
        start_index: Optional[int] = 0,
        end_index: Optional[int] = 1
    ):
    filename, ext = os.path.splitext(os.path.basename(filepath))
    if ext in ACCEPT_TIFF_EXT:
        filepath = join_filepath([DIRPATH.INPUT_DIR, filepath])
        save_dirpath = join_filepath([
            os.path.dirname(filepath),
            filename,
        ])
        json_filepath = join_filepath([
            save_dirpath,
            f'{filename}_{str(start_index)}_{str(end_index)}.json'
        ])
        if not os.path.exists(json_filepath):
            save_tiff2json(filepath, save_dirpath, start_index, end_index)
    else:
        json_filepath = filepath

    return JsonReader.read_as_output(json_filepath)

@router.post("/outputs/image/{filepath:path}/add_roi", response_model=EditRoiSuccess)
async def add_roi(filepath: str, pos: RoiPos):
    if not os.path.exists(filepath):
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND)
    if os.path.basename(filepath) != 'cell_roi.json':
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST)
    
    try:
        cell_roi_data, max_index = execute_add_ROI(node_dirpath=os.path.dirname(filepath), pos=[pos.posx, pos.posy, pos.sizex, pos.sizey])
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))
    return EditRoiSuccess(data=cell_roi_data, max_index=max_index)

@router.post("/outputs/image/{filepath:path}/merge_roi", response_model=EditRoiSuccess)
async def merge_roi(filepath: str, roi_list: RoiList):
    if not os.path.exists(filepath):
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND)
    if os.path.basename(filepath) != 'cell_roi.json':
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST)
    
    try:
        cell_roi_data, max_index = execute_merge_roi(node_dirpath=os.path.dirname(filepath), merged_roi_ids=roi_list.ids)
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))
    return EditRoiSuccess(data=cell_roi_data, max_index=max_index)
    
@router.post("/outputs/image/{filepath:path}/delete_roi", response_model=EditRoiSuccess)
async def delete_roi(filepath: str, roi_list: RoiList):
    if not os.path.exists(filepath):
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND)
    if os.path.basename(filepath) != 'cell_roi.json':
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST)
    
    try:
        cell_roi_data, max_index = excute_delete_roi(node_dirpath=os.path.dirname(filepath), delete_roi_ids=roi_list.ids)
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))
    return EditRoiSuccess(data=cell_roi_data, max_index=max_index)


@router.get("/outputs/csv/{filepath:path}")
async def get_csv(filepath: str):
    filepath = join_filepath([DIRPATH.INPUT_DIR, filepath])

    filename, _ = os.path.splitext(os.path.basename(filepath))
    save_dirpath = join_filepath([
        os.path.dirname(filepath),
        filename
    ])
    create_directory(save_dirpath)
    json_filepath = join_filepath([
        save_dirpath,
        f'{filename}.json'
    ])

    JsonWriter.write_as_split(
        json_filepath,
        pd.read_csv(filepath, header=None)
    )
    return JsonReader.read_as_output(json_filepath)
