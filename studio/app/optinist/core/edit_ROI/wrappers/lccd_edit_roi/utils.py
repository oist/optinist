from studio.app.optinist.core.nwb.nwb import NWBDATASET


def set_nwbfile(lccd_data, roi_list, function_id, fluorescence=None):
    nwbfile = {}

    nwbfile[NWBDATASET.ROI] = {function_id: roi_list}

    if fluorescence is not None:
        nwbfile[NWBDATASET.FLUORESCENCE] = {
            function_id: {
                "Fluorescence": {
                    "table_name": "ROIs",
                    "region": list(range(len(fluorescence))),
                    "name": "Fluorescence",
                    "data": fluorescence,
                    "unit": "lumens",
                }
            }
        }

    nwbfile[NWBDATASET.COLUMN] = {
        function_id: {
            "name": "iscell",
            "description": "two columns - iscell & probcell",
            "data": lccd_data.get("is_cell"),
        }
    }

    # NWB追加
    nwbfile[NWBDATASET.POSTPROCESS] = {
        function_id: {
            "add_roi": lccd_data.get("add_roi", []),
            "delete_roi": lccd_data.get("delete_roi", []),
            "merge_roi": lccd_data.get("merge_roi", []),
        }
    }

    return nwbfile
