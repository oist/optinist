from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import join_filepath
from optinist.routers.model import FILETYPE


def smk_input(config, name):
    output_list = []
    for value in config["rules"].values():
        if value["type"] == name:
            if value["type"] in [FILETYPE.IMAGE]:
                output_list.append([
                    join_filepath([DIRPATH.INPUT_DIR, x])
                    for x in value["input"]
                ])
            elif value["type"] in [FILETYPE.CSV, FILETYPE.BEHAVIOR, FILETYPE.HDF5]:
                output_list.append(join_filepath([DIRPATH.INPUT_DIR, value["input"]]))
            else:
                output_list.append([
                    join_filepath([
                        DIRPATH.OUTPUT_DIR, x
                    ]) for x in value["input"]
                ])
    return output_list


def smk_output(config, name):
    output_list = []
    for key, value in config["rules"].items():
        if value["type"] == name:
            output_list.append(join_filepath([DIRPATH.OUTPUT_DIR, value["output"]]))

    return output_list
