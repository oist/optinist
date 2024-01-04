from typing import List

import h5py
import numpy as np
from fastapi import APIRouter

from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.dir_path import DIRPATH
from studio.app.optinist.schemas.hdf5 import HDF5Node

router = APIRouter()


class HDF5Getter:
    @classmethod
    def get(cls, filepath) -> List[HDF5Node]:
        cls.hdf5_list = []
        with h5py.File(filepath, "r") as f:
            f.visititems(cls.get_ds_dictionaries)

        return cls.hdf5_list

    @classmethod
    def get_ds_dictionaries(cls, path: str, node: h5py.Dataset):
        if isinstance(node, h5py.Dataset):
            if len(node.shape) != 0:
                cls.recursive_dir_tree(cls.hdf5_list, path.split("/"), node, "")

    @classmethod
    def recursive_dir_tree(
        cls,
        node_list: List[HDF5Node],
        path_list: List[str],
        node: h5py.Dataset,
        parent_path: str,
    ):
        name = path_list[0]
        if name.startswith("#"):
            return

        path = name if parent_path == "" else f"{parent_path}/{name}"

        is_exists = False
        # 既にkeyがある
        for i, value in enumerate(node_list):
            if value.name == name:
                is_exists = True
                if len(path_list) > 1:
                    cls.recursive_dir_tree(
                        node_list[i].nodes, path_list[1:], node, path
                    )

        if not is_exists:
            if len(path_list) > 1:
                node_list.append(
                    HDF5Node(
                        isDir=True,
                        name=name,
                        path=path,
                        nodes=[],
                    )
                )
                cls.recursive_dir_tree(node_list[-1].nodes, path_list[1:], node, path)
            else:
                node_list.append(
                    HDF5Node(
                        isDir=False,
                        name=name,
                        path=path,
                        shape=node.shape,
                        nbytes=f"{int(node.nbytes / (1000**2))} M",
                        dataType="array"
                        if isinstance(node[:], np.ndarray)
                        else type(node[:]).__name__,
                    )
                )


@router.get("/hdf5/{file_path:path}", response_model=List[HDF5Node], tags=["outputs"])
async def get_files(file_path: str, workspace_id: str):
    file_path = join_filepath([DIRPATH.INPUT_DIR, workspace_id, file_path])
    return HDF5Getter.get(file_path)
