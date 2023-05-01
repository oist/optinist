import os

from optinist.api.dataclass.dataclass import ImageData, Suite2pData
from optinist.api.utils.filepath_creater import create_directory, join_filepath


def suite2p_file_convert(
    image: ImageData, output_dir: str, params: dict = None
) -> dict(ops=Suite2pData):
    from suite2p import default_ops, io

    print("start suite2_file_convert")

    data_path_list = []
    data_name_list = []
    for file_path in image.path:
        data_path_list.append(os.path.dirname(file_path))
        data_name_list.append(os.path.basename(file_path))

    print(data_path_list)
    print(data_name_list)
    # data pathと保存pathを指定
    db = {
        "data_path": data_path_list,
        "tiff_list": data_name_list,
        "save_path0": output_dir,
        "save_folder": "suite2p",
    }

    ops = {**default_ops(), **params, **db}

    # save folderを指定
    create_directory(join_filepath([ops["save_path0"], ops["save_folder"]]))

    # save ops.npy(parameter) and data.bin
    ops = io.tiff_to_binary(ops.copy())

    info = {
        "meanImg": ImageData(
            ops["meanImg"], output_dir=output_dir, file_name="meanImg"
        ),
        "ops": Suite2pData(ops, file_name="ops"),
    }

    return info
