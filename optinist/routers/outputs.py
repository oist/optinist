from typing import Optional
from fastapi import APIRouter
from glob import glob
import os
import json
from dataclasses import dataclass
from typing import Dict

from optinist.cui_api.json_writer import save_tiff2json, save_csv2json
from optinist.cui_api.dir_path import DIRPATH
from optinist.cui_api.filepath_creater import join_filepath

router = APIRouter()


class JsonReader:
    @classmethod
    def json_read(cls, filepath):
        with open(filepath, 'r') as f:
            json_data = json.load(f)
        return json_data

    @classmethod
    def timeseries_read(cls, filepath):
        json_data = cls.json_read(filepath)
        return JsonTimeSeriesData(
            xrange=list(json_data["data"].keys()),
            data=json_data["data"],
            std=json_data["std"] if "std" in json_data else None,
        )


@dataclass
class JsonData:
    data: Dict[str, dict]


@dataclass
class JsonTimeSeriesData(JsonData):
    xrange: list
    # data: Dict[str, dict]
    std: Dict[str, dict]


@router.get("/outputs/timedata/{dirpath:path}")
async def read_file(dirpath: str, index: Optional[int] = None):
    json_data = JsonReader.timeseries_read(join_filepath([dirpath, f'{str(index)}.json']))

    return_data = JsonTimeSeriesData(
        xrange=[],
        data={},
        std={},
    )

    if index == 0:
        num_files = len(glob(join_filepath([dirpath, '*.json'])))
        data = {
            str(i): {
                json_data.xrange[0]: json_data.data[json_data.xrange[0]]
            }
            for i in range(num_files)
        }

        std = {}
        if json_data.std is not None:
            std = {
                str(i): {
                    json_data.xrange[0]: json_data.data[json_data.xrange[0]]
                }
                for i in range(num_files)
            }

        return_data = JsonTimeSeriesData(
            xrange=json_data.xrange,
            data=data,
            std=std,
        )

    str_index = str(index)
    return_data.data[str_index] = json_data.data
    if json_data.std is not None:
        return_data.std[str_index] = json_data.std

    return return_data


@router.get("/outputs/alltimedata/{dirpath:path}")
async def read_file(dirpath: str):
    return_data = JsonTimeSeriesData(
        xrange=[],
        data={},
        std={},
    )
    for i, path in enumerate(glob(join_filepath([dirpath, '*.json']))):
        str_i = str(i)
        json_data = JsonReader.timeseries_read(path)
        if i == 0:
            return_data.xrange = json_data.xrange

        return_data.data[str_i] = json_data.data
        if json_data.std is not None:
            return_data.std[str_i] = json_data.std

    return return_data


@router.get("/outputs/data/{filepath:path}")
async def read_file(filepath: str):
    return JsonData(JsonReader.json_read(filepath))


@router.get("/outputs/html/{file_path:path}")
async def read_file(file_path: str):
    with open(file_path, 'r') as f:
        data = f.read()

    return { "data": data }


@router.get("/outputs/image/{file_path:path}")
async def read_image(
        file_path: str,
        start_index: Optional[int] = None,
        end_index: Optional[int] = None
    ):
    file_name, ext = os.path.splitext(os.path.basename(file_path))
    if ext in ['.tif', '.tiff', '.TIF', '.TIFF']:
        folder_path = os.path.dirname(file_path)
        tiff_file_path = file_path
        file_path = join_filepath(
            [folder_path, f'{file_name}_{str(start_index)}_{str(end_index)}.json'])
        if not os.path.exists(file_path):
            save_tiff2json(tiff_file_path, start_index, end_index)
    elif ext == '.json':
        pass
    else:
        assert False, "Extension Error"

    with open(file_path, 'r') as f:
        json_dict = json.load(f)

    return { "data": json_dict }


@router.get("/outputs/csv/{file_path:path}")
async def read_csv(file_path: str):
    file_name, ext = os.path.splitext(os.path.basename(file_path))

    if ext == '.csv':
        folder_path = os.path.dirname(file_path)
        csv_file_path = file_path
        file_path = join_filepath([folder_path, f'{file_name}.json'])

        save_csv2json(csv_file_path)

    with open(file_path, 'r') as f:
        json_dict = json.load(f)

    return { "data": json_dict }
