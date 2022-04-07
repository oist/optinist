import os
import shutil
from glob import glob
import sys

if sys.version_info >= (3, 8):
    from typing import List, Optional, TypedDict
else:
    from typing import List, Optional
    from typing_extensions import TypedDict

from fastapi import APIRouter, File, Response, UploadFile, Form
from optinist.cui_api.const import BASE_DIR
from optinist.cui_api.utils import join_file_path

router = APIRouter()

ACCEPT_FILE_TYPES = ["tif", "tiff", "TIF", "TIFF", "json", "csv", "hdf5", "nwb"]


class TreeNode(TypedDict):
    path: str
    name: str
    isdir: bool
    nodes: Optional[List["TreeNode"]]


def get_accept_files(path: str, file_types: List[str]):
    files_list = []
    for file_type in file_types:
        files_list.extend(glob(
            join_file_path([path, "**", f"*.{file_type}"]), recursive=True))

    return files_list


def get_dir_tree(dir_path: str, file_types: List[str]) -> List[TreeNode]:
    nodes: List[TreeNode] = []
    for node_name in os.listdir(dir_path):
        node_path = join_file_path([dir_path, node_name])
        if os.path.isfile(node_path) and node_name.endswith(tuple(file_types)):
            nodes.append({
                "path": node_path,
                "name": node_name,
                "isdir": False,
            })
        elif os.path.isdir(node_path) and len(get_accept_files(node_path, file_types)) > 0:
            nodes.append({
                "path": node_path,
                "name": node_name,
                "isdir": True,
                "nodes": get_dir_tree(node_path, file_types)
            })

    return nodes


@router.get("/files")
async def get_files(file_type: Optional[str] = None):
    tree = []
    if(file_type is None):
        tree = get_dir_tree(BASE_DIR, ACCEPT_FILE_TYPES)
    else:
        if file_type == "image":
            tree = get_dir_tree(BASE_DIR, ["tif", 'TIF', 'tiff', 'TIFF'])
        elif file_type == "csv":
            tree = get_dir_tree(BASE_DIR, ["csv"])
        elif file_type == "hdf5":
            tree = get_dir_tree(BASE_DIR, ["hdf5", "nwb"])
        else:
            # TODO 他のファイル種別の仕様が分かり次第追加
            pass

    return tree


@router.post("/files/upload/{fileName}")
async def create_file(response: Response, fileName: str, file: UploadFile = File(...)):
    root_dir = os.path.splitext(join_file_path([BASE_DIR, fileName]))[0]
    if not os.path.exists(root_dir):
        os.makedirs(root_dir, exist_ok=True)

    file_path = join_file_path([root_dir, fileName])

    with open(file_path, "wb") as f:
        shutil.copyfileobj(file.file, f)

    return { "file_path": file_path }
