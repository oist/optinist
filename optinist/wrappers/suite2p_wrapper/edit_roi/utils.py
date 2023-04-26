import time

import numpy as np
from scipy import stats

from optinist.api.dataclass.dataclass import *
from optinist.api.nwb.nwb import NWBDATASET
from optinist.api.nwb.nwb_creater import overwrite_nwb


def masks_and_traces(ops, stat_manual, stat_orig):
    """main extraction function
    inputs: ops and stat
    creates cell and neuropil masks and extracts traces
    returns: F (ROIs x time), Fneu (ROIs x time), F_chan2, Fneu_chan2, ops, stat
    F_chan2 and Fneu_chan2 will be empty if no second channel
    """
    from suite2p import detection, extraction

    if "aspect" in ops:
        dy, dx = int(ops["aspect"] * 10), 10
    else:
        d0 = ops["diameter"]
        dy, dx = (d0, d0) if isinstance(d0, int) else d0
    t0 = time.time()
    # Concatenate stat so a good neuropil function can be formed
    stat_all = stat_orig
    for n in range(len(stat_manual)):
        stat_all.append(stat_manual[n])
    stat_all = detection.stats.roi_stats(stat_all, dy, dx, ops["Ly"], ops["Lx"])

    cell_masks = [
        extraction.masks.create_cell_mask(
            stat, Ly=ops["Ly"], Lx=ops["Lx"], allow_overlap=ops["allow_overlap"]
        )
        for stat in stat_all
    ]
    cell_pix = extraction.masks.create_cell_pix(stat_all, Ly=ops["Ly"], Lx=ops["Lx"])
    manual_roi_stats = stat_all[-len(stat_manual) :]
    manual_cell_masks = cell_masks[-len(stat_manual) :]

    manual_neuropil_masks = extraction.masks.create_neuropil_masks(
        ypixs=[stat["ypix"] for stat in manual_roi_stats],
        xpixs=[stat["xpix"] for stat in manual_roi_stats],
        cell_pix=cell_pix,
        inner_neuropil_radius=ops["inner_neuropil_radius"],
        min_neuropil_pixels=ops["min_neuropil_pixels"],
    )
    print("Masks made in %0.2f sec." % (time.time() - t0))

    F, Fneu, F_chan2, Fneu_chan2, ops = extraction.extract_traces_from_masks(
        ops, manual_cell_masks, manual_neuropil_masks
    )

    # compute activity statistics for classifier
    npix = np.array([stat_orig[n]["npix"] for n in range(len(stat_orig))]).astype(
        "float32"
    )
    for n in range(len(manual_roi_stats)):
        manual_roi_stats[n]["npix_norm"] = manual_roi_stats[n]["npix"] / np.mean(
            npix[:100]
        )  # What if there are less than 100 cells?
        manual_roi_stats[n]["compact"] = 1
        manual_roi_stats[n]["footprint"] = 2
        manual_roi_stats[n]["manual"] = 1  # Add manual key

    # subtract neuropil and compute skew, std from F
    dF = F - ops["neucoeff"] * Fneu
    sk = stats.skew(dF, axis=1)
    sd = np.std(dF, axis=1)

    for n in range(F.shape[0]):
        manual_roi_stats[n]["skew"] = sk[n]
        manual_roi_stats[n]["std"] = sd[n]
        manual_roi_stats[n]["med"] = [
            np.mean(manual_roi_stats[n]["ypix"]),
            np.mean(manual_roi_stats[n]["xpix"]),
        ]

    dF = extraction.preprocess(
        F=dF,
        baseline=ops["baseline"],
        win_baseline=ops["win_baseline"],
        sig_baseline=ops["sig_baseline"],
        fs=ops["fs"],
        prctile_baseline=ops["prctile_baseline"],
    )

    spks = extraction.dcnv.oasis(
        F=dF, batch_size=ops["batch_size"], tau=ops["tau"], fs=ops["fs"]
    )

    return F, Fneu, F_chan2, Fneu_chan2, spks, ops, manual_roi_stats


def get_stat0_add_roi(ops, pos):
    posx, posy, sizex, sizey = list(map(round, pos))
    xrange = (np.arange(-1 * int(sizex), 1) + int(posx)).astype(np.int32)
    yrange = (np.arange(-1 * int(sizey), 1) + int(posy)).astype(np.int32)
    xrange += int(np.floor(sizex / 2))
    yrange += int(np.floor(sizey / 2))

    ellipse = np.zeros((yrange.size, xrange.size), np.bool)
    x, y = np.meshgrid(np.arange(0, xrange.size, 1), np.arange(0, yrange.size, 1))
    ellipse = (
        (y - sizey / 2) ** 2 / (sizey / 2) ** 2
        + (x - sizex / 2) ** 2 / (sizex / 2) ** 2
    ) <= 1

    ellipse = ellipse[:, np.logical_and(xrange >= 0, xrange < ops["Lx"])]
    xrange = xrange[np.logical_and(xrange >= 0, xrange < ops["Lx"])]
    ellipse = ellipse[np.logical_and(yrange >= 0, yrange < ops["Ly"]), :]
    yrange = yrange[np.logical_and(yrange >= 0, yrange < ops["Ly"])]

    med = pos[2:4]
    x, y = np.meshgrid(xrange, yrange)
    ypix = y[ellipse].flatten()
    xpix = x[ellipse].flatten()
    lam = np.ones(ypix.shape)

    return [{"ypix": ypix, "xpix": xpix, "lam": lam, "npix": ypix.size, "med": med}]


def save_json_data(ops, im, save_path, save_data=[]):
    info = {}
    iscell = ops.get("iscell")

    for d in save_data:
        if d == "ops":
            info["ops"] = Suite2pData(ops)

        if d == "max_proj":
            info["max_proj"] = ImageData(ops["max_proj"], file_name="max_proj")

        if d == "Vcorr":
            info["Vcorr"] = ImageData(ops["Vcorr"], file_name="Vcorr")

        if d == "fluorescence":
            info["fluorescence"] = FluoData(ops["F"], file_name="fluorescence")

        if d == "iscell" and iscell is not None:
            info["iscell"] = IscellData(ops["iscell"], file_name="iscell")

        if d == "all_roi" and im is not None:
            info["all_roi"] = RoiData(np.nanmax(im, axis=0), file_name="all_roi")

        if d == "non_cell_roi" and iscell is not None and im is not None:
            info["non_cell_roi"] = RoiData(
                np.nanmax(im[~iscell], axis=0), file_name="noncell_roi"
            )

        if d == "cell_roi" and iscell is not None and im is not None:
            info["cell_roi"] = RoiData(
                np.nanmax(im[iscell], axis=0), file_name="cell_roi"
            )

        if d == "nwbfile":
            nwbfile = set_nwbfile(ops)
            overwrite_nwb(nwbfile, save_path, "suite2p_roi.nwb")

    for k, v in info.items():
        if isinstance(v, BaseData):
            v.save_json(save_path)

    return info


def set_nwbfile(ops):
    stat = ops.get("stat")
    iscell = ops.get("iscell")
    F = ops.get("F")
    Fneu = ops.get("Fneu")

    roi_list = []
    for i in range(len(stat)):
        kargs = {}
        kargs["pixel_mask"] = np.array(
            [stat[i]["ypix"], stat[i]["xpix"], stat[i]["lam"]]
        ).T
        roi_list.append(kargs)
    nwbfile = {}

    nwbfile[NWBDATASET.ROI] = {"roi_list": roi_list}

    # iscellを追加
    nwbfile[NWBDATASET.COLUMN] = {
        "roi_column": {
            "name": "iscell",
            "discription": "two columns - iscell & probcell",
            "data": iscell,
        }
    }

    # Fluorenceを追加
    nwbfile[NWBDATASET.FLUORESCENCE] = {}
    for name, data in zip(["Fluorescence", "Neuropil"], [F, Fneu]):
        nwbfile[NWBDATASET.FLUORESCENCE][name] = {
            "table_name": name,
            "region": list(range(len(data))),
            "name": name,
            "data": data,
            "unit": "lumens",
            "rate": ops["fs"],
        }

    add_roi = ops.get("add_roi") if ops.get("add_roi") else []
    delete_roi = ops.get("delete_roi") if ops.get("delete_roi") else []
    merge_roi = ops.get("merge_roi") if ops.get("merge_roi") else []

    # NWB追加
    nwbfile[NWBDATASET.POSTPROCESS] = {
        "add_roi": add_roi,
        "delete_roi": delete_roi,
        "merge_roi": merge_roi,
    }

    return nwbfile
