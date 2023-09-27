import os
import shutil
from glob import glob
from typing import List

from fastapi import APIRouter, Depends, File, UploadFile

from studio.app.common.core.utils.filepath_creater import (
    create_directory,
    join_filepath,
)
from studio.app.common.core.workspace.workspace_dependencies import (
    is_workspace_available,
    is_workspace_owner,
)
from studio.app.common.schemas.files import FilePath, TreeNode
from studio.app.const import ACCEPT_CSV_EXT, ACCEPT_HDF5_EXT, ACCEPT_TIFF_EXT, FILETYPE
from studio.app.dir_path import DIRPATH

router = APIRouter(prefix="/files", tags=["files"])


class DirTreeGetter:
    @classmethod
    def get_tree(
        cls, workspace_id, file_types: List[str], dirname: str = None
    ) -> List[TreeNode]:
        nodes: List[TreeNode] = []

        if dirname is None:
            absolute_dirpath = join_filepath([DIRPATH.INPUT_DIR, workspace_id])
        else:
            absolute_dirpath = join_filepath([DIRPATH.INPUT_DIR, workspace_id, dirname])

        if not os.path.exists(absolute_dirpath):
            return nodes

        sorted_listdir = sorted(
            os.listdir(absolute_dirpath),
            key=lambda x: (not os.path.isdir(join_filepath([absolute_dirpath, x])), x),
        )
        for node_name in sorted_listdir:
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
                        nodes=cls.get_tree(workspace_id, file_types, relative_path),
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


@router.get(
    "/{workspace_id}",
    response_model=List[TreeNode],
    dependencies=[Depends(is_workspace_available)],
)
async def get_files(workspace_id: str, file_type: str = None):
    if file_type == FILETYPE.IMAGE:
        return DirTreeGetter.get_tree(workspace_id, ACCEPT_TIFF_EXT)
    elif file_type == FILETYPE.CSV:
        return DirTreeGetter.get_tree(workspace_id, ACCEPT_CSV_EXT)
    elif file_type == FILETYPE.HDF5:
        return DirTreeGetter.get_tree(workspace_id, ACCEPT_HDF5_EXT)
    else:
        return []


@router.post(
    "/{workspace_id}/upload/{filename}",
    response_model=FilePath,
    dependencies=[Depends(is_workspace_owner)],
)
async def create_file(workspace_id: str, filename: str, file: UploadFile = File(...)):
    create_directory(join_filepath([DIRPATH.INPUT_DIR, workspace_id]))

    filepath = join_filepath([DIRPATH.INPUT_DIR, workspace_id, filename])

    with open(filepath, "wb") as f:
        shutil.copyfileobj(file.file, f)

    return {"file_path": filename}
