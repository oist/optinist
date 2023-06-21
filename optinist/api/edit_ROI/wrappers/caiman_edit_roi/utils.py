from optinist.api.nwb.nwb import NWBDATASET


def set_nwbfile(cnmf_data):
    im = cnmf_data.get("im")
    is_cell = cnmf_data.get("is_cell")
    fluorescence = cnmf_data.get("fluorescence")

    # NWBの追加
    nwbfile = {}
    # NWBにROIを追加
    roi_list = []
    n_cells = im.shape[0]
    for i in range(n_cells):
        kargs = {}
        kargs["image_mask"] = im[i, :]
        roi_list.append(kargs)

    nwbfile[NWBDATASET.ROI] = {"roi_list": roi_list}

    # iscellを追加
    nwbfile[NWBDATASET.COLUMN] = {
        "roi_column": {
            "name": "iscell",
            "discription": "two columns - iscell & probcell",
            "data": is_cell,
        }
    }

    # Fluorescence

    nwbfile[NWBDATASET.FLUORESCENCE] = {
        "RoiResponseSeries": {
            "table_name": "ROIs",
            "region": list(range(n_cells)),
            "name": "RoiResponseSeries",
            "data": fluorescence.T,
            "unit": "lumens",
        },
    }

    # NWB追加
    nwbfile[NWBDATASET.POSTPROCESS] = {
        "add_roi": cnmf_data.get("add_roi", []),
        "delete_roi": cnmf_data.get("delete_roi", []),
        "merge_roi": cnmf_data.get("merge_roi", []),
    }

    return nwbfile
