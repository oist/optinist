from typing import Optional
from fastapi import APIRouter
from glob import glob
import os
import json
from dataclasses import dataclass
from typing import Dict

from optinist.api.utils.json_writer import save_tiff2json, save_csv2json
from optinist.api.utils.filepath_creater import join_filepath

router = APIRouter()


@dataclass
class Data:
    data: Dict[str, dict]


@dataclass
class JsonTimeSeriesData(Data):
    xrange: list
    std: Dict[str, dict]


class Reader:
    @classmethod
    def read(cls, filepath):
        with open(filepath, 'r') as f:
            data = f.read()
        return Data(data=data)


class JsonReader:
    @classmethod
    def read(cls, filepath):
        with open(filepath, 'r') as f:
            json_data = json.load(f)
        return json_data

    @classmethod
    def timeseries_read(cls, filepath) -> JsonTimeSeriesData:
        json_data = cls.read(filepath)
        return JsonTimeSeriesData(
            xrange=list(json_data["data"].keys()),
            data=json_data["data"],
            std=json_data["std"] if "std" in json_data else None,
        )


@router.get("/outputs/timedata/{dirpath:path}")
async def read_file(dirpath: str, index: int):
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
            std=std if json_data.std is not None else {},
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
    return Data(JsonReader.read(filepath))


@router.get("/outputs/html/{filepath:path}")
async def read_html_file(filepath: str):
    return Data(Reader.read(filepath))


@router.get("/outputs/image/{filepath:path}")
async def read_image(
        filepath: str,
        start_index: Optional[int] = None,
        end_index: Optional[int] = None
    ):
    filename, ext = os.path.splitext(os.path.basename(filepath))
    if ext in ['.tif', '.tiff', '.TIF', '.TIFF']:
        json_filepath = join_filepath(
            [os.path.dirname(filepath), f'{filename}_{str(start_index)}_{str(end_index)}.json'])
        if not os.path.exists(json_filepath):
            save_tiff2json(filepath, start_index, end_index)
    else:
        json_filepath = filepath

    return Data(JsonReader.read(json_filepath))


@router.get("/outputs/csv/{filepath:path}")
async def read_csv(filepath: str):
    dirpath = os.path.dirname(filepath)
    filename, _ = os.path.splitext(os.path.basename(filepath))
    json_filepath = join_filepath([dirpath, f'{filename}.json'])
    save_csv2json(filepath, json_filepath)
    return Data(JsonReader.read(json_filepath))
