from optinist.api.dataclass.dataclass import CaimanCnmfData, FluoData, RoiData
from optinist.api.edit_ROI.utils import create_mask
from optinist.api.edit_ROI.wrappers.caiman_edit_roi.utils import set_nwbfile


def execute_add_ROI(node_dirpath, posx, posy, sizex, sizey):
    import numpy as np

    # load data
    cnmf_data = np.load(f"{node_dirpath}/caiman_cnmf.npy", allow_pickle=True).item()
    is_cell = cnmf_data.get("is_cell")
    im = cnmf_data.get("im")
    images = cnmf_data.get("images")
    add_roi = cnmf_data.get("add_roi", [])
    fluorescence = cnmf_data.get("fluorescence")

    # Create mask for new roi
    new_roi = create_mask(posx, posy, sizex, sizey, images.shape[1:])

    im = np.concatenate((im, new_roi[np.newaxis, :, :]), axis=0)
    is_cell = np.append(is_cell, True)

    # extract fluorescence of new ROI
    num_frames = images.shape[0]
    reshapeImages = images.reshape([images.shape[2] * images.shape[1], num_frames])
    new_fluorescence = np.zeros(num_frames)
    new_fluorescence = np.mean(reshapeImages[new_roi.reshape(-1, 1)[:, 0] > 0], axis=0)
    fluorescence = np.vstack([fluorescence, new_fluorescence])

    cell_roi = np.zeros(im.shape)
    num_rois = im.shape[0]
    for i in range(num_rois):
        cell_roi[i, :, :] = np.where(im[i, :, :] != 0, i + 1, np.nan)
    add_roi.append(num_rois)

    # save data
    cnmf_data["im"] = im
    cnmf_data["is_cell"] = is_cell
    cnmf_data["fluorescence"] = fluorescence
    cnmf_data["add_roi"] = add_roi

    info = {
        "fluorescence": FluoData(fluorescence, file_name="fluorescence"),
        "cell_roi": RoiData(
            np.nanmax(cell_roi[is_cell], axis=0),
            output_dir=node_dirpath,
            file_name="cell_roi",
        ),
        "cnmf_data": CaimanCnmfData(cnmf_data),
        "nwbfile": set_nwbfile(cnmf_data),
    }

    return info
