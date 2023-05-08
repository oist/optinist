import os

import numpy as np

from optinist.api.dataclass.dataclass import FluoData, RoiData, Suite2pData
from optinist.api.edit_ROI.wrappers.suite2p_edit_roi.utils import set_nwbfile


def execute_merge_roi(node_dirpath, ids):
    from suite2p import ROI
    from suite2p.detection.stats import median_pix, roi_stats

    ops = np.load(os.path.join(node_dirpath, "suite2p.npy"), allow_pickle=True).item()
    iscell = ops.get("iscell")
    stat_orig = ops.get("stat")
    im = ops.get("im")
    F_all = ops["F"]
    Fneu_all = ops["Fneu"]
    merge_roi = ops.get("merge_roi", [])

    print("merging activity... this may take some time")
    ypix = np.zeros((0,), np.int32)
    xpix = np.zeros((0,), np.int32)
    lam = np.zeros((0,), np.float32)
    footprints = np.array([])
    F = np.zeros((0, F_all.shape[1]), np.float32)
    Fneu = np.zeros((0, F_all.shape[1]), np.float32)

    # merged cells id
    merged_cells = []
    for n in np.array(ids):
        merged_cells.append(n)
    merged_cells = np.unique(np.array(merged_cells))

    # merge cells regions
    for n in merged_cells:
        # deactivate merged cells
        iscell[n] = False
        ypix = np.append(ypix, stat_orig[n]["ypix"])
        xpix = np.append(xpix, stat_orig[n]["xpix"])
        lam = np.append(lam, stat_orig[n]["lam"])
        footprints = np.append(footprints, stat_orig[n]["footprint"])
        F = np.append(F, F_all[n, :][np.newaxis, :], axis=0)
        Fneu = np.append(Fneu, Fneu_all[n, :][np.newaxis, :], axis=0)

    # remove overlaps from merged cells regions
    ipix = np.concatenate((ypix[:, np.newaxis], xpix[:, np.newaxis]), axis=1)
    _, goodi = np.unique(ipix, return_index=True, axis=0)
    ypix = ypix[goodi]
    xpix = xpix[goodi]
    lam = lam[goodi]

    # compute statistics of merges
    stat0 = {}
    stat0["ypix"] = ypix
    stat0["xpix"] = xpix
    stat0["lam"] = lam / lam.sum()
    stat0["npix"] = ypix.size
    stat0["med"] = median_pix(ypix, xpix)

    stat0["chan2_prob"] = -1
    stat0["inmerge"] = -1

    # compute activity of merged cells
    F = F.mean(axis=0)
    Fneu = Fneu.mean(axis=0)

    if "aspect" in ops:
        d0 = np.array([int(ops["aspect"] * 10), 10])
    else:
        d0 = ops["diameter"]
        if isinstance(d0, int):
            d0 = [d0, d0]

    # add cell to structs
    print(type(stat_orig), type(stat0), type(stat_orig[0]))
    stat_orig.append(stat0)
    stat_orig = roi_stats(stat_orig, d0[0], d0[1], ops["Ly"], ops["Lx"])
    stat_orig[-1]["lam"] = stat_orig[-1]["lam"] * merged_cells.size

    array = np.expand_dims(
        ROI(
            ypix=stat_orig[-1]["ypix"],
            xpix=stat_orig[-1]["xpix"],
            lam=stat_orig[-1]["lam"],
            med=stat_orig[-1]["med"],
            do_crop=False,
        ).to_array(Ly=ops["Ly"], Lx=ops["Lx"]),
        axis=0,
    )
    id = np.nanmax(im) + 1
    array *= id
    merge_roi.append(float(id))
    merge_roi += (merged_cells + 1).tolist()
    merge_roi.append((-1.0))

    im = np.concatenate((im, array), axis=0)
    im[im == 0] = np.nan

    try:
        F_all = np.concatenate((ops["F"], np.expand_dims(F, axis=0)), axis=0)
        Fneu_all = np.concatenate((ops["Fneu"], np.expand_dims(Fneu, axis=0)), axis=0)
        iscell = np.append(iscell, True)
    except Exception as e:
        print(e)

    ops["stat"] = stat_orig
    ops["F"] = F_all
    ops["Fneu"] = Fneu_all
    ops["iscell"] = iscell
    ops["im"] = im
    ops["merge_roi"] = merge_roi

    info = {
        "ops": Suite2pData(ops),
        "fluorescence": FluoData(ops["F"], file_name="fluorescence"),
        "all_roi": RoiData(
            np.nanmax(im, axis=0), output_dir=node_dirpath, file_name="all_roi"
        ),
        "non_cell_roi": RoiData(
            np.nanmax(im[~iscell], axis=0),
            output_dir=node_dirpath,
            file_name="noncell_roi",
        ),
        "cell_roi": RoiData(
            np.nanmax(im[iscell], axis=0), output_dir=node_dirpath, file_name="cell_roi"
        ),
        "nwbfile": set_nwbfile(ops),
    }
    return info
