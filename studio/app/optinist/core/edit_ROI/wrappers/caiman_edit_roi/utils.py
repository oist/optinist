from studio.app.optinist.core.nwb.nwb import NWBDATASET
from studio.app.optinist.dataclass.roi import EditRoiData


def set_nwbfile(edit_roi_data: EditRoiData, iscell, function_id, fluorescence=None):
    # NWBの追加
    nwbfile = {}

    # NWBにROIを追加
    roi_list = []
    n_cells = edit_roi_data.im.shape[0]
    for i in range(n_cells):
        kargs = {}
        kargs["image_mask"] = edit_roi_data.im[i, :]
        roi_list.append(kargs)
    nwbfile[NWBDATASET.ROI] = {function_id: roi_list}

    nwbfile[NWBDATASET.COLUMN] = {
        function_id: {
            "name": "iscell",
            "description": "two columns - iscell & probcell",
            "data": iscell,
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
            "add_roi": edit_roi_data.add_roi,
            "delete_roi": edit_roi_data.delete_roi,
            "merge_roi": edit_roi_data.merge_roi,
        }
    }

    return nwbfile
