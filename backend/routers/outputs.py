from fastapi import APIRouter
import os
import json

router = APIRouter()

@router.get("/outputs/{file_path:path}")
async def read_file(file_path: str):
    with open(os.path.join(".", file_path), 'r') as f:
        json_dict = json.load(f)
    return { "data": json_dict }
