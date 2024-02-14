import numpy as np

from studio.app.common.dataclass.image import ImageData
from studio.app.optinist.dataclass.iscell import IscellData
from studio.app.optinist.dataclass.roi import RoiData


def roi_from_hdf5(
    cell_img: ImageData,
    output_dir: str,
    iscell: IscellData = None,
    params: dict = None,
    **kwargs
) -> dict():
    if iscell is not None:
        iscell = iscell.data
        data = (
            cell_img.data[iscell != 0]
            if params.get("show_iscell", True)
            else cell_img.data[iscell == 0]
        )
        return {
            "roi": RoiData(
                np.nanmax(data, axis=0),
                output_dir=output_dir,
                file_name="roi",
            ),
            "nwbfile": {},
        }
    else:
        return {
            "roi": RoiData(
                np.nanmax(cell_img.data, axis=0),
                output_dir=output_dir,
                file_name="roi",
            ),
            "nwbfile": {},
        }
