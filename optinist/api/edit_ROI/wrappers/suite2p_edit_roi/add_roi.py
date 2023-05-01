import os

import numpy as np

from optinist.api.dataclass.dataclass import FluoData, RoiData, Suite2pData
from optinist.api.edit_ROI.wrappers.suite2p_edit_roi.utils import (
    get_stat0_add_roi,
    masks_and_traces,
    set_nwbfile,
)


def execute_add_ROI(node_dirpath, posx, posy, sizex, sizey):
    from suite2p import ROI

    ops = np.load(os.path.join(node_dirpath, "suite2p.npy"), allow_pickle=True).item()
    iscell = ops.get("iscell")
    stat = ops.get("stat")
    stat0 = get_stat0_add_roi(ops, posx, posy, sizex, sizey)
    im = ops.get("im")
    add_roi = ops.get("add_roi", [])

    Fcell, Fneu, _, _, _, ops, manual_stat_cell = masks_and_traces(ops, stat0, stat)

    for n in range(len(manual_stat_cell)):
        array = np.expand_dims(
            ROI(
                ypix=manual_stat_cell[n]["ypix"],
                xpix=manual_stat_cell[n]["xpix"],
                lam=manual_stat_cell[n]["lam"],
                med=manual_stat_cell[n]["med"],
                do_crop=False,
            ).to_array(Ly=ops["Ly"], Lx=ops["Lx"]),
            axis=0,
        )
        id = np.nanmax(im) + 1
        array *= id
        add_roi.append(id)
        im = np.concatenate((im, array), axis=0)
    im[im == 0] = np.nan

    # Save fluorescence traces
    F_all = np.concatenate((ops["F"], Fcell), axis=0)
    Fneu_all = np.concatenate((ops["Fneu"], Fneu), axis=0)
    iscell = np.append(iscell, True)
    ops["stat"] = stat
    ops["F"] = F_all
    ops["Fneu"] = Fneu_all
    ops["iscell"] = iscell
    ops["im"] = im
    ops["add_roi"] = add_roi

    info = {
        "ops": Suite2pData(ops),
        "fluorescence": FluoData(ops["F"], file_name="fluorescence"),
        "all_roi": RoiData(
            np.nanmax(im, axis=0), output_dir=node_dirpath, file_name="all_roi"
        ),
        "cell_roi": RoiData(
            np.nanmax(im[iscell], axis=0), output_dir=node_dirpath, file_name="cell_roi"
        ),
        "nwbfile": set_nwbfile(ops),
    }
    return info
