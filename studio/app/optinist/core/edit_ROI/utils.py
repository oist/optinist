import os
from concurrent.futures import ProcessPoolExecutor
from typing import Tuple

import numpy as np
from fastapi import HTTPException, status
from snakemake import snakemake

from studio.app.common.core.logger import LOGGING_CLIENT_ID_KEY, AppLogger
from studio.app.common.core.logger_context_helpers import (
    get_client_id_for_subprocess,
    with_client_id_context,
)
from studio.app.common.core.utils.filepath_finder import find_condaenv_filepath
from studio.app.dir_path import DIRPATH
from studio.app.optinist.core.edit_ROI.wrappers import edit_roi_wrapper_dict
from studio.app.optinist.schemas.roi import RoiPos

logger = AppLogger.get_logger()


class EditRoiUtils:
    @classmethod
    def conda(cls, config):
        algo = config["algo"]
        if "conda_name" in edit_roi_wrapper_dict[algo]:
            conda_name = edit_roi_wrapper_dict[algo]["conda_name"]
            return find_condaenv_filepath(conda_name) if conda_name else None

        return None

    @classmethod
    def get_algo(cls, filepath):
        algo_list = edit_roi_wrapper_dict.keys()

        algo = next((algo for algo in algo_list if algo in filepath), None)
        if not algo:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST)
        return algo

    @classmethod
    def execute(cls, filepath: set):
        client_id = get_client_id_for_subprocess()

        result = False

        with ProcessPoolExecutor(max_workers=1) as executor:
            logger.info("start snakemake edit_roi process.")

            future = executor.submit(
                cls._execute_process, filepath, client_id=client_id
            )
            result = future.result()

            logger.info("finish snakemake edit_roi process. result: %s", result)

        if not result:
            logger.error("edit_ROI snakemake run failed.")
            raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR)

    @classmethod
    @with_client_id_context  # Automatically set client_id for logging
    def _execute_process(cls, filepath: str, client_id: str = None) -> bool:
        result = snakemake(
            DIRPATH.SNAKEMAKE_EDIT_ROI_FILEPATH,
            use_conda=True,
            cores=2,
            workdir=f"{os.path.dirname(DIRPATH.STUDIO_DIR)}",
            config={
                "type": "EDIT_ROI",
                "algo": cls.get_algo(filepath),
                "file_path": filepath,
                LOGGING_CLIENT_ID_KEY: client_id,
            },
        )

        return result


def create_ellipse_mask(shape: Tuple[int, int], roi_pos: RoiPos):
    x, y, width, height = (
        round(roi_pos.posx),
        round(roi_pos.posy),
        round(roi_pos.sizex),
        round(roi_pos.sizey),
    )

    x_coords = np.arange(0, shape[0])
    y_coords = np.arange(0, shape[1])
    xx, yy = np.meshgrid(x_coords, y_coords)

    # Calculate the distance of each pixel from the center of the ellipse
    a = width / 2
    b = height / 2
    distance = ((xx - x) / a) ** 2 + ((yy - y) / b) ** 2

    # Set the pixels within the ellipse to 1 and the pixels outside to NaN
    ellipse = np.empty(shape)
    ellipse[:] = np.nan
    ellipse[distance <= 1] = 1

    return ellipse
