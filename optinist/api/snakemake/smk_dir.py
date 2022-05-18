from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import join_filepath
from optinist.routers.model import FILETYPE


def smk_input(config, name):
    for value in config["rules"].values():
        if value["type"] == name:
            if value["type"] in [FILETYPE.IMAGE]:
                return [
                    join_filepath([DIRPATH.INPUT_DIR, x])
                    for x in value["input"]
                ]
            elif value["type"] in [FILETYPE.CSV, FILETYPE.BEHAVIOR, FILETYPE.HDF5]:
                return join_filepath([DIRPATH.INPUT_DIR, value["input"]])
            else:
                return [
                    join_filepath([DIRPATH.OUTPUT_DIR, x])
                    for x in value["input"]
                ]


def smk_output(config, name):
    for value in config["rules"].values():
        if value["type"] == name:
            return join_filepath([DIRPATH.OUTPUT_DIR, value["output"]])
