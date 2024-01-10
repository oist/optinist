from functools import reduce
from typing import List

import numpy as np
from fastapi import APIRouter
from pymatreader import read_mat

from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.dir_path import DIRPATH
from studio.app.optinist.schemas.mat import MatNode

router = APIRouter()


class MatGetter:
    @classmethod
    def data(cls, filepath, dataPath: str = None):
        data = read_mat(filepath)
        data = {
            key: value
            for key, value in data.items()
            if not key.startswith("__") and not key.startswith("#")
        }
        keys = dataPath.split("/") if dataPath is not None else []
        return reduce(lambda d, key: d[key], keys, data)

    @classmethod
    def get(cls, filepath, workspace_id) -> List[MatNode]:
        filepath = join_filepath([DIRPATH.INPUT_DIR, workspace_id, filepath])
        data = cls.data(filepath, dataPath=None)
        return [cls.dict_to_matnode(value, key, key) for key, value in data.items()]

    @classmethod
    def dict_to_matnode(cls, data, name, current_path=""):
        if isinstance(data, dict):
            return MatNode(
                isDir=True,
                name=name,
                path=current_path,
                nodes=[
                    cls.dict_to_matnode(
                        v, k, f"{current_path}/{k}" if current_path else k
                    )
                    for k, v in data.items()
                ],
            )
        elif isinstance(data, np.ndarray):
            return MatNode(
                isDir=False,
                name=name,
                path=current_path,
                shape=data.shape,
                dataType="array",
                nbytes=f"{int(data.nbytes / (1000**2))} M",
            )
        else:
            return MatNode(
                isDir=False,
                name=name,
                path=current_path,
                dataType=type(data).__name__,
            )


@router.get(
    "/mat/{file_path:path}",
    response_model=List[MatNode],
    tags=["outputs"],
)
async def get_matfiles(file_path: str, workspace_id: str):
    return MatGetter.get(file_path, workspace_id)
