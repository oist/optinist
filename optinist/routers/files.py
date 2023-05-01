import os
import shutil
from glob import glob
from typing import List

from fastapi import APIRouter, File, UploadFile

from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import create_directory, join_filepath
from optinist.routers.const import ACCEPT_CSV_EXT, ACCEPT_HDF5_EXT, ACCEPT_TIFF_EXT
from optinist.routers.model import FILETYPE, FilePath, TreeNode

router = APIRouter()


class DirTreeGetter:
    @classmethod
    def get_tree(cls, file_types: List[str], dirname: str = None) -> List[TreeNode]:
        nodes: List[TreeNode] = []

        if dirname is None:
            absolute_dirpath = DIRPATH.INPUT_DIR
        else:
            absolute_dirpath = join_filepath([DIRPATH.INPUT_DIR, dirname])

        for node_name in os.listdir(absolute_dirpath):
            if dirname is None:
                relative_path = node_name
            else:
                relative_path = join_filepath([dirname, node_name])

            search_dirpath = join_filepath([absolute_dirpath, node_name])

            if os.path.isfile(search_dirpath) and node_name.endswith(tuple(file_types)):
                nodes.append(
                    TreeNode(
                        path=relative_path,
                        name=node_name,
                        isdir=False,
                        nodes=[],
                    )
                )
            elif (
                os.path.isdir(search_dirpath)
                and len(cls.accept_files(search_dirpath, file_types)) > 0
            ):
                nodes.append(
                    TreeNode(
                        path=node_name,
                        name=node_name,
                        isdir=True,
                        nodes=cls.get_tree(file_types, relative_path),
                    )
                )

        return nodes

    @classmethod
    def accept_files(cls, path: str, file_types: List[str]):
        files_list = []
        for file_type in file_types:
            files_list.extend(
                glob(join_filepath([path, "**", f"*{file_type}"]), recursive=True)
            )

        return files_list


@router.get("/files", response_model=List[TreeNode], tags=["files"])
async def get_files(file_type: str = None):
    if file_type == FILETYPE.IMAGE:
        return DirTreeGetter.get_tree(ACCEPT_TIFF_EXT)
    elif file_type == FILETYPE.CSV:
        return DirTreeGetter.get_tree(ACCEPT_CSV_EXT)
    elif file_type == FILETYPE.HDF5:
        return DirTreeGetter.get_tree(ACCEPT_HDF5_EXT)


@router.post("/files/upload/{filename}", response_model=FilePath, tags=["files"])
async def create_file(filename: str, file: UploadFile = File(...)):
    create_directory(DIRPATH.INPUT_DIR)

    filepath = join_filepath([DIRPATH.INPUT_DIR, filename])

    with open(filepath, "wb") as f:
        shutil.copyfileobj(file.file, f)

    return {"file_path": filename}
