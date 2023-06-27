import os
from glob import glob
from typing import Optional

from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import join_filepath


def find_filepath(name, category) -> Optional[str]:
    name, _ = os.path.splitext(name)

    filepaths = glob(
        join_filepath([DIRPATH.ROOT_DIR, "wrappers", "**", category, f"{name}.yaml"]),
        recursive=True,
    )
    return filepaths[0] if len(filepaths) > 0 else None


def find_param_filepath(name: str):
    if name in ["snakemake", "nwb"]:
        return join_filepath([DIRPATH.ROOT_DIR, "api", name, f"{name}.yaml"])
    else:
        return find_filepath(name, "params")


def find_condaenv_filepath(name: str):
    return find_filepath(name, "env")
