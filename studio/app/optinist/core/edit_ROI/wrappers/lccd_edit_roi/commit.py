import numpy as np

from studio.app.common.dataclass.image import ImageData
from studio.app.optinist.core.edit_ROI.wrappers.lccd_edit_roi.utils import set_nwbfile
from studio.app.optinist.dataclass import EditRoiData, FluoData, IscellData, RoiData


def commit_edit(
    images: ImageData,
    data: EditRoiData,
    fluorescence: FluoData,
    iscell,
    node_dirpath,
    function_id,
):
    from studio.app.optinist.core.edit_ROI.edit_ROI import CellType

    fluorescence = fluorescence.data
    num_cell = data.im.shape[0]

    new_fluorescences = np.zeros((num_cell, fluorescence.shape[1]))
    new_fluorescences[: len(fluorescence)] = fluorescence

    iscell[iscell == CellType.TEMP_DELETE] = CellType.NON_ROI
    for i in range(num_cell):
        if iscell[i] == CellType.TEMP_ADD:
            new_fluorescences[i] = np.mean(images[:, ~np.isnan(data.im[i])], axis=1)
            iscell[i] = CellType.ROI

    data.commit()

    info = {
        "cell_roi": RoiData(
            np.nanmax(data.im[iscell != CellType.NON_ROI], axis=0),
            output_dir=node_dirpath,
            file_name="cell_roi",
        ),
        "fluorescence": FluoData(new_fluorescences, file_name="fluorescence"),
        "iscell": IscellData(iscell),
        "edit_roi_data": data,
        "nwbfile": set_nwbfile(data, iscell, function_id, new_fluorescences),
    }
    return info
