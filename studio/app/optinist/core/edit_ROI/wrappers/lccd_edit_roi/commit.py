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

    # Special processing if saving with no ROIs
    if all(i == CellType.TEMP_DELETE for i in iscell):
        empty_fluorescences = np.zeros((0, fluorescence.shape[1]))  # Keep time dim
        iscell[iscell == CellType.TEMP_DELETE] = CellType.NON_ROI
        data.commit()

        info = {
            "cell_roi": RoiData(
                np.full(data.im.shape[1:], np.nan),
                output_dir=node_dirpath,
                file_name="cell_roi",
            ),
            "fluorescence": FluoData(empty_fluorescences, file_name="fluorescence"),
            "iscell": IscellData(iscell),
            "edit_roi_data": data,
            "nwbfile": set_nwbfile(data, iscell, function_id, empty_fluorescences),
        }
        return info

    if fluorescence.size > 0:
        new_fluorescences = np.zeros((num_cell, fluorescence.shape[1]))
        new_fluorescences[: len(fluorescence)] = fluorescence
    else:
        # If starting with no ROIs, initialize with correct time dimension
        new_fluorescences = np.zeros((num_cell, images.shape[0]))

    num_committed = len(fluorescence)
    added = iscell == CellType.TEMP_ADD
    iscell[iscell == CellType.TEMP_DELETE] = CellType.NON_ROI
    for i in range(num_committed, num_cell):
        new_fluorescences[i] = np.mean(images[:, ~np.isnan(data.im[i])], axis=1)
        if added[i]:
            iscell[i] = CellType.ROI

    data.commit()

    info = {
        "cell_roi": RoiData(
            np.nanmax(data.im[iscell != CellType.NON_ROI], axis=0, initial=np.nan),
            output_dir=node_dirpath,
            file_name="cell_roi",
        ),
        "fluorescence": FluoData(new_fluorescences, file_name="fluorescence"),
        "iscell": IscellData(iscell),
        "edit_roi_data": data,
        "nwbfile": set_nwbfile(data, iscell, function_id, new_fluorescences),
    }
    return info
