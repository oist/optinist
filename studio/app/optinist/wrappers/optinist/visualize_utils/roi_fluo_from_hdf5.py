import numpy as np

from studio.app.common.core.logger import AppLogger
from studio.app.common.dataclass.image import ImageData
from studio.app.optinist.dataclass.fluo import FluoData
from studio.app.optinist.dataclass.iscell import IscellData
from studio.app.optinist.dataclass.roi import RoiData
from studio.app.optinist.wrappers.optinist.visualize_utils.roi_from_hdf5 import (
    max_projection,
)


def roi_fluo_from_hdf5(
    cell_img: ImageData,
    fluo: FluoData,
    output_dir: str,
    iscell: IscellData = None,
    params: dict = None,
    **kwargs,
) -> dict():
    """
    Processes ROI and fluorescence data, aligning them, and adds enriched metadata.

    Parameters:
        cell_img (ImageData): Cell image data containing ROI information.
        fluo (FluoData): Fluorescence data.
        output_dir (str): Directory to save the output data.
        iscell (IscellData): Binary data indicating cell regions (optional).
        params (dict): Optional parameters for additional processing.
        **kwargs: Additional keyword arguments.

    Returns:
        dict: A dictionary containing processed ROI and fluorescence data.
    """
    logger = AppLogger.get_logger()
    logger.info("Starting ROI and fluorescence data processing")

    # Base output data
    all_roi = RoiData(
        np.nanmax(cell_img.data, axis=0), output_dir=output_dir, file_name="all_roi"
    )

    if params["transpose"]:
        fluorescence = FluoData(np.transpose(fluo.data), file_name="fluorescence")
    else:
        fluorescence = FluoData(fluo.data, file_name="fluorescence")

    if iscell is None:
        return {"all_roi": all_roi, "fluorescence": fluorescence}

    # Extract iscell data
    iscell_data = iscell.data

    # Process ROI types
    non_cell_roi = RoiData(
        max_projection(cell_img.data[np.where(iscell_data == 0)]),
        output_dir=output_dir,
        file_name="noncell_roi",
    )

    cell_roi = RoiData(
        max_projection(cell_img.data[np.where(iscell_data != 0)]),
        output_dir=output_dir,
        file_name="cell_roi",
    )

    # Construct output dictionary
    output = {
        "iscell": IscellData(iscell_data, file_name="iscell"),
        "all_roi": all_roi,
        "non_cell_roi": non_cell_roi,
        "cell_roi": cell_roi,
        "fluorescence": fluorescence,
    }

    logger.info("Completed ROI and fluorescence data processing")
    return output
