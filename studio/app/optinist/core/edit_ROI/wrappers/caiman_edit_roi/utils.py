from studio.app.optinist.core.nwb.nwb import NWBDATASET


def set_nwbfile(cnmf_data, function_id):
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

    nwbfile[NWBDATASET.ROI] = {function_id: roi_list}

    # iscellを追加
    nwbfile[NWBDATASET.COLUMN] = {
        function_id: {
            "name": "iscell",
            "discription": "two columns - iscell & probcell",
            "data": is_cell,
        }
    }

    # Fluorescence

    nwbfile[NWBDATASET.FLUORESCENCE] = {
        function_id: {
            "Fluorescence": {
                "table_name": "ROIs",
                "region": list(range(n_cells)),
                "name": "Fluorescence",
                "data": fluorescence.T,
                "unit": "lumens",
            },
        }
    }

    # NWB追加
    nwbfile[NWBDATASET.POSTPROCESS] = {
        function_id: {
            "add_roi": cnmf_data.get("add_roi", []),
            "delete_roi": cnmf_data.get("delete_roi", []),
            "merge_roi": cnmf_data.get("merge_roi", []),
        }
    }

    return nwbfile
