from typing import Optional
from fastapi import APIRouter
import os
import json
from .utils.save import save_tiff2json, save_csv2json
from .const import BASE_DIR

router = APIRouter()


@router.get("/outputs/data/{file_path:path}")
async def read_file(file_path: str):
    with open(os.path.join(BASE_DIR, file_path), 'r') as f:
        json_dict = json.load(f)
    return { "data": json_dict }


@router.get("/outputs/image/{file_path:path}")
async def read_image(file_path: str, max_index: Optional[int] = None):
    file_path = os.path.join(BASE_DIR, file_path)
    file_name, ext = os.path.splitext(os.path.basename(file_path))

    if ext == '.tif' or ext == '.TIF':
        folder_path = os.path.dirname(file_path)
        tiff_file_path = file_path
        file_path = os.path.join(
            folder_path, f'{file_name}_{str(max_index)}.json')

        if not os.path.exists(file_path):
            save_tiff2json(tiff_file_path, max_index)

    with open(file_path, 'r') as f:
        json_dict = json.load(f)

    return { "data": json_dict}


@router.get("/outputs/csv/{file_path:path}")
async def read_csv(file_path: str):
    file_path = os.path.join(BASE_DIR, file_path)
    file_name, ext = os.path.splitext(os.path.basename(file_path))

    if ext == '.csv':
        folder_path = os.path.dirname(file_path)
        csv_file_path = file_path
        file_path = os.path.join(folder_path, f'{file_name}.json')

        if not os.path.exists(file_path):
            save_csv2json(csv_file_path)

    with open(file_path, 'r') as f:
        json_dict = json.load(f)
    return { "data": json_dict}
