import numpy as np

from studio.app.optinist.core.nwb.nwb import NWBDATASET
from studio.app.optinist.dataclass.roi import EditRoiData


def set_nwbfile(output_info: dict, function_id):
    edit_roi_data: EditRoiData = output_info.get("edit_roi_data")
    iscell = output_info.get("iscell").data

    # NWBの追加
    nwbfile = {}

    # NWBにROIを追加
    nwbfile[NWBDATASET.ROI] = {function_id: _get_roilist(output_info, function_id)}

    # breakpoint()
    nwbfile[NWBDATASET.COLUMN] = {
        function_id: {
            "name": "iscell",
            "discription": "two columns - iscell & probcell",
            "data": iscell,
        }
    }

    # Fluorescence
    nwbfile[NWBDATASET.FLUORESCENCE] = _get_fluor(output_info, function_id)

    # NWB追加
    nwbfile[NWBDATASET.POSTPROCESS] = {
        function_id: {
            "add_roi": edit_roi_data.add_roi,
            "delete_roi": edit_roi_data.delete_roi,
            "merge_roi": edit_roi_data.merge_roi,
        }
    }

    return nwbfile


def _get_roilist(output_info: dict, function_id):
    roi_list = []
    # if "suite2p" in function_id:
    #     breakpoint()
    #     stat = output_info.get("ops").data['stat']
    #     for i in stat:
    #         kargs = {}
    #         kargs["pixel_mask"] = np.array([i["ypix"], i["xpix"], i["lam"]]).T
    #         roi_list.append(kargs)
    # else:
    im = output_info.get("edit_roi_data").im
    n_cells = im.shape[0]
    for i in range(n_cells):
        kargs = {}
        kargs["image_mask"] = im[i, :]
        roi_list.append(kargs)
    return roi_list


def _get_fluor(output_info: dict, function_id):
    # Fluorenceを追加
    if "suite2p" in function_id:
        ops = output_info.get("ops").data
        F = ops["F"]
        Fneu = ops["Fneu"]
        return {
            function_id: {
                "Fluorescence": {
                    "table_name": "Fluorescence",
                    "region": list(range(len(F))),
                    "name": "Fluorescence",
                    "data": F,
                    "unit": "lumens",
                    "rate": ops["fs"],
                },
                "Neuropil": {
                    "table_name": "Neuropil",
                    "region": list(range(len(Fneu))),
                    "name": "Neuropil",
                    "data": Fneu,
                    "unit": "lumens",
                    "rate": ops["fs"],
                },
            }
        }
    else:
        fluorescence = output_info.get("fluorescence").data
        return {
            function_id: {
                "Fluorescence": {
                    "table_name": "ROIs",
                    "region": list(range(len(fluorescence))),
                    "name": "Fluorescence",
                    "data": fluorescence,
                    "unit": "lumens",
                },
            }
        }
