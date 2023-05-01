import os

import numpy as np

from optinist.wrappers.suite2p_wrapper.edit_roi.utils import save_json_data


def excute_delete_roi(node_dirpath, delete_roi_ids):
    ops = np.load(os.path.join(node_dirpath, "suite2p.npy"), allow_pickle=True).item()
    iscell = ops.get("iscell")
    im = ops.get("im")

    delete_roi = ops.get("delete_roi") if ops.get("delete_roi") else []
    for id in delete_roi_ids:
        iscell[id] = False
        delete_roi.append(id + 1)

    ops["iscell"] = iscell
    ops["delete_roi"] = delete_roi

    save_json_data(
        ops,
        im,
        save_path=node_dirpath,
        save_data=["ops", "non_cell_roi", "cell_roi", "nwbfile"],
    )

    max_index = len(ops["F"])
    return max_index
