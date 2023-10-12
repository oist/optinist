from typing import Tuple

import numpy as np

from studio.app.optinist.core.nwb.nwb import NWBDATASET
from studio.app.optinist.dataclass.roi import EditRoiData
from studio.app.optinist.schemas.roi import RoiPos


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


def create_ellipse_mask(shape: Tuple[int, int], roi_pos: RoiPos):
    x, y, width, height = (
        round(roi_pos.posx),
        round(roi_pos.posy),
        round(roi_pos.sizex),
        round(roi_pos.sizey),
    )

    x_coords = np.arange(0, shape[0])
    y_coords = np.arange(0, shape[1])
    xx, yy = np.meshgrid(x_coords, y_coords)

    # Calculate the distance of each pixel from the center of the ellipse
    a = width / 2
    b = height / 2
    distance = ((xx - x) / a) ** 2 + ((yy - y) / b) ** 2

    # Set the pixels within the ellipse to 1 and the pixels outside to NaN
    ellipse = np.empty(shape)
    ellipse[:] = np.nan
    ellipse[distance <= 1] = 1

    return ellipse
