import os

from studio.app.common.core.experiment.experiment import ExptOutputPathIds
from studio.app.common.core.logger import AppLogger
from studio.app.common.core.utils.filepath_creater import (
    create_directory,
    join_filepath,
)
from studio.app.common.dataclass import ImageData
from studio.app.optinist.dataclass import Suite2pData

logger = AppLogger.get_logger()


def suite2p_file_convert(
    image: ImageData, output_dir: str, params: dict = None, **kwargs
) -> dict(ops=Suite2pData):
    import numpy as np
    import tifffile
    from suite2p import default_ops, io
    from suite2p.io.tiff import open_tiff, use_sktiff_reader

    function_id = ExptOutputPathIds(output_dir).function_id
    logger.info("start suite2p_file_convert: %s", function_id)

    data_path_list = []
    data_name_list = []
    for file_path in image.path:
        use_sktiff = (
            True
            if params.get("force_sktiff")
            else use_sktiff_reader(file_path, batch_size=params.get("batch_size"))
        )

        tif, _ = open_tiff(file_path, use_sktiff)
        if use_sktiff:
            im = tifffile.imread(file_path)
        else:
            im = tif.data()

        if im.dtype.type == np.float32:
            im = (im - im.min()) / (im.max() - im.min())
            im = (im * np.iinfo(np.int16).max).astype(np.int16)

        file_path = join_filepath([output_dir, os.path.basename(file_path)])
        tifffile.imwrite(file_path, im)

        data_path_list.append(os.path.dirname(file_path))
        data_name_list.append(os.path.basename(file_path))

    logger.info(data_path_list)
    logger.info(data_name_list)
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
