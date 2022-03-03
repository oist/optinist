from typing import Optional, List
from fastapi import APIRouter
import os
import json
from glob import glob
from .utils.save import save_tiff2json, save_csv2json
from cui_api.const import BASE_DIR
from cui_api.utils import join_file_path

router = APIRouter()


@router.get("/outputs/timedata/{file_path:path}")
async def read_file(file_path: str, index: Optional[int] = None):
    _dir = file_path
    
    with open(join_file_path([_dir, f'{str(index)}.json']), 'r') as f:
        json_dict = json.load(f)

    return_dict = {}
    if index == 0:
        num_files = len(glob(join_file_path([_dir, '*.json'])))
        return_dict = {str(i): {0: json_dict["0"]["0"]} for i in range(num_files)}

    return_dict[str(index)] = json_dict["0"]

    return { "data": return_dict }


@router.get("/outputs/alltimedata/{file_path:path}")
async def read_file(file_path: str):
    _dir = file_path

    return_dict = {}
    for index, path in enumerate(glob(join_file_path([_dir, '*.json']))):
        with open(path, 'r') as f:
            json_dict = json.load(f)
            return_dict[str(index)] = json_dict["0"]

    return { "data": return_dict }


@router.get("/outputs/data/{file_path:path}")
async def read_file(file_path: str):
    with open(file_path, 'r') as f:
        json_dict = json.load(f)

    return { "data": json_dict }


@router.get("/outputs/image/{file_path:path}")
async def read_image(
        file_path: str,
        start_index: Optional[int] = None,
        end_index: Optional[int] = None
    ):
    # file_path = join_file_path([BASE_DIR, file_path])
    file_name, ext = os.path.splitext(os.path.basename(file_path))
    if ext in ['.tif', '.tiff', '.TIF', '.TIFF']:
        folder_path = os.path.dirname(file_path)
        tiff_file_path = file_path
        file_path = join_file_path(
            [folder_path, f'{file_name}_{str(start_index)}_{str(end_index)}.json'])
        if not os.path.exists(file_path):
            save_tiff2json(tiff_file_path, start_index, end_index)
    elif ext == '.json':
        pass
    else:
        assert False, "Extension Error"

    with open(file_path, 'r') as f:
        json_dict = json.load(f)

    return { "data": json_dict}


@router.get("/outputs/csv/{file_path:path}")
async def read_csv(file_path: str):
    # file_path = join_file_path([BASE_DIR, file_path])
    file_name, ext = os.path.splitext(os.path.basename(file_path))

    if ext == '.csv':
        folder_path = os.path.dirname(file_path)
        csv_file_path = file_path
        file_path = join_file_path([folder_path, f'{file_name}.json'])

        # if not os.path.exists(file_path):
        save_csv2json(csv_file_path)

    with open(file_path, 'r') as f:
        json_dict = json.load(f)

    return { "data": json_dict}
