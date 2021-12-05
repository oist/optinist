from fastapi import APIRouter
import os
import json
from .utils import save_tiff_to_json

router = APIRouter()

@router.get("/outputs/{file_path:path}")
async def read_file(file_path: str):
    file_path = os.path.join(".", file_path)
    print(file_path)

    file_name, ext = os.path.splitext(os.path.basename(file_path))

    if ext == '.tif':
        folder_path = os.path.dirname(file_path)
        file_path = os.path.join(folder_path, f'{file_name}.json')
        if not os.path.exists(file_path):
            save_tiff_to_json(file_path, 10)

    with open(file_path, 'r') as f:
        json_dict = json.load(f)

    return { "data": json_dict}
