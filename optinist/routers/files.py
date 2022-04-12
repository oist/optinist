import os
import shutil
from glob import glob
from typing import List
from fastapi import APIRouter, File, UploadFile

from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import join_filepath
from optinist.routers.model import FILETYPE, TreeNode

router = APIRouter()


ACCEPT_TIFF_EXT = ["tif", "tiff", "TIF", "TIFF"]
ACCEPT_CSV_EXT = ["csv"]
ACCEPT_HDF5_EXT = ["hdf5", "nwb", "HDF5", "NWB"]


def get_accept_files(path: str, file_types: List[str]):
    files_list = []
    for file_type in file_types:
        files_list.extend(glob(
            join_filepath([path, "**", f"*.{file_type}"]), recursive=True))

    return files_list


def get_dir_tree(dirpath: str, file_types: List[str]) -> List[TreeNode]:
    nodes: List[TreeNode] = []
    for node_name in os.listdir(dirpath):
        node_path = join_filepath([dirpath, node_name])
        if os.path.isfile(node_path) and node_name.endswith(tuple(file_types)):
            nodes.append(TreeNode(
                path=node_path,
                name=node_name,
                isdir=False,
                nodes=[],
            ))
        elif os.path.isdir(node_path) and len(get_accept_files(node_path, file_types)) > 0:
            nodes.append(TreeNode(
                path=node_path,
                name=node_name,
                isdir=True,
                nodes=get_dir_tree(node_path, file_types)
            ))

    return nodes


@router.get("/files")
async def get_files(file_type: str = None):
    if file_type == FILETYPE.IMAGE:
        return get_dir_tree(DIRPATH.BASE_DIR, ACCEPT_TIFF_EXT)
    elif file_type == FILETYPE.CSV:
        return get_dir_tree(DIRPATH.BASE_DIR, ACCEPT_CSV_EXT)
    elif file_type == FILETYPE.HDF5:
        return get_dir_tree(DIRPATH.BASE_DIR, ACCEPT_HDF5_EXT)


@router.post("/files/upload/{filename}")
async def create_file(filename: str, file: UploadFile = File(...)):
    dirpath = os.path.splitext(join_filepath([DIRPATH.BASE_DIR, filename]))[0]
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)

    filepath = join_filepath([dirpath, filename])

    with open(filepath, "wb") as f:
        shutil.copyfileobj(file.file, f)

    return { "file_path": filepath }
