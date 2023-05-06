from optinist.api.nwb.nwb import NWBDATASET


def set_nwbfile(lccd_data, roi_list, fluorescence=None):
    nwbfile = {}

    nwbfile[NWBDATASET.ROI] = {"roi_list": roi_list}

    if fluorescence is not None:
        nwbfile[NWBDATASET.FLUORESCENCE] = {}
        nwbfile[NWBDATASET.FLUORESCENCE]["Fluorescence"] = {
            "table_name": "Fluorescence",
            "region": list(range(len(fluorescence))),
            "name": "Fluorescence",
            "data": fluorescence,
            "unit": "lumens",
        }

    nwbfile[NWBDATASET.COLUMN] = {
        "roi_column": {
            "name": "iscell",
            "discription": "two columns - iscell & probcell",
            "data": lccd_data.get("is_cell"),
        }
    }

    # NWB追加
    nwbfile[NWBDATASET.POSTPROCESS] = {
        "add_roi": lccd_data.get("add_roi", []),
        "delete_roi": lccd_data.get("delete_roi", []),
        "merge_roi": lccd_data.get("merge_roi", []),
    }

    return nwbfile
