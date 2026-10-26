import numpy as np

from studio.app.common.dataclass.image import ImageData
from studio.app.optinist.dataclass.iscell import IscellData
from studio.app.optinist.dataclass.roi import RoiData


def max_projection(rois):
    # not nanmax(initial=np.nan): it raises on the integer stacks an HDF5 can hold
    return np.nanmax(rois, axis=0) if len(rois) else np.full(rois.shape[1:], np.nan)


def roi_from_hdf5(
    cell_img: ImageData,
    output_dir: str,
    iscell: IscellData = None,
    params: dict = None,
    **kwargs
) -> dict(iscell=IscellData):
    if iscell is not None:
        iscell_data = iscell.data
        return {
            "iscell": IscellData(iscell_data, file_name="iscell"),
            "all_roi": RoiData(
                np.nanmax(cell_img.data, axis=0),
                output_dir=output_dir,
                file_name="roi",
            ),
            "non_cell_roi": RoiData(
                max_projection(cell_img.data[iscell_data == 0]),
                output_dir=output_dir,
                file_name="noncell_roi",
            ),
            "cell_roi": RoiData(
                max_projection(cell_img.data[iscell_data != 0]),
                output_dir=output_dir,
                file_name="cell_roi",
            ),
        }
    else:
        return {
            "all_roi": RoiData(
                np.nanmax(cell_img.data, axis=0),
                output_dir=output_dir,
                file_name="roi",
            ),
        }
