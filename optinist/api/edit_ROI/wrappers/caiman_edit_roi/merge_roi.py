from optinist.api.dataclass.dataclass import CaimanCnmfData, FluoData, RoiData
from optinist.api.edit_ROI.wrappers.caiman_edit_roi.utils import set_nwbfile


def execute_merge_roi(node_dirpath: str, ids: list):
    import numpy as np

    # load data
    cnmf_data = np.load(f"{node_dirpath}/caiman_cnmf.npy", allow_pickle=True).item()
    is_cell = cnmf_data.get("is_cell")
    merge_roi = cnmf_data.get("merge_roi", [])
    im = cnmf_data.get("im")
    fluorescence = cnmf_data.get("fluorescence")

    # get merging ROI
    merging_ROIs = []
    [merging_ROIs.append(im[id, :, :]) for id in ids]

    # get merged ROI
    merged_ROI = np.maximum.reduce(merging_ROIs)
    is_cell[ids] = False
    is_cell = np.append(is_cell, True)
    im = np.concatenate((im, merged_ROI[np.newaxis, :, :]), axis=0)

    # get merged F
    merged_f = np.mean((fluorescence[ids, :]), axis=0)
    fluorescence = np.vstack([fluorescence, merged_f])

    cell_roi = np.zeros(im.shape)
    num_rois = im.shape[0]
    for i in range(num_rois):
        cell_roi[i, :, :] = np.where(im[i, :, :] != 0, i + 1, np.nan)
    merge_roi.append(float(num_rois))
    merge_roi += [(id + 1) for id in ids]
    merge_roi.append((-1.0))

    cnmf_data["im"] = im
    cnmf_data["is_cell"] = is_cell
    cnmf_data["fluorescence"] = fluorescence
    cnmf_data["merge_roi"] = merge_roi

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
