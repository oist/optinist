import os

from optinist.api.dataclass.dataclass import FluoData, LccdData, RoiData
from optinist.api.edit_ROI.utils import create_mask
from optinist.api.edit_ROI.wrappers.lccd_edit_roi.utils import set_nwbfile


def execute_add_ROI(node_dirpath, posx, posy, sizex, sizey):
    import numpy as np

    lccd_data = np.load(
        os.path.join(node_dirpath, "lccd.npy"), allow_pickle=True
    ).item()

    images = lccd_data.get("images", None)
    roi = lccd_data.get("roi", None)
    is_cell = lccd_data.get("is_cell")
    add_roi = lccd_data.get("add_roi", [])

    shape = images.shape[:2]
    num_frames = images.shape[2]
    dff_f0_frames = 100
    dff_f0_percentile = 8

    new_roi = create_mask(posx, posy, sizex, sizey, shape).reshape(-1, 1)
    roi = np.hstack((roi, new_roi))
    is_cell = np.append(is_cell, True)

    num_cell = roi.shape[1]
    images = images.reshape([images.shape[0] * images.shape[1], images.shape[2]])
    timeseries = np.zeros([num_cell, num_frames])

    add_roi.append(num_cell)

    # Get ROI list and extract Fluor
    im = []
    for i in range(num_cell):
        im.append((roi[:, i].reshape(shape)) * (i + 1))
        timeseries[i, :] = np.mean(images[roi[:, i] > 0, :], axis=0)
    im = np.stack(im)
    im[im == 0] = np.nan

    # Get DFF
    timeseries_dff = np.ones([num_cell, num_frames]) * np.nan
    for i in range(num_cell):
        for k in range(num_frames):
            if (k - dff_f0_frames >= 0) and (k + dff_f0_frames < num_frames):
                f0 = np.percentile(
                    timeseries[i, k - dff_f0_frames : k + dff_f0_frames],
                    dff_f0_percentile,
                )
                timeseries_dff[i, k] = (timeseries[i, k] - f0) / f0

    roi_list = [{"image_mask": roi[:, i].reshape(shape)} for i in range(num_cell)]

    lccd_data["roi"] = roi
    lccd_data["is_cell"] = is_cell
    lccd_data["add_roi"] = add_roi

    info = {
        "lccd": LccdData(lccd_data),
        "cell_roi": RoiData(
            np.nanmax(im[is_cell], axis=0),
            output_dir=node_dirpath,
            file_name="cell_roi",
        ),
        "fluorescence": FluoData(timeseries, file_name="fluorescence"),
        "dff": FluoData(timeseries_dff, file_name="dff"),
        "nwbfile": set_nwbfile(lccd_data, roi_list, fluorescence=timeseries),
    }

    return info
