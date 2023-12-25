import numpy as np

from studio.app.optinist.core.edit_ROI.wrappers.suite2p_edit_roi.utils import (
    get_stat0_add_roi,
    masks_and_traces,
    set_nwbfile,
)
from studio.app.optinist.dataclass import (
    EditRoiData,
    FluoData,
    IscellData,
    RoiData,
    Suite2pData,
)
from studio.app.optinist.schemas.roi import RoiPos


def commit_edit(data: EditRoiData, ops: Suite2pData, iscell, node_dirpath, function_id):
    from suite2p.detection.stats import median_pix, roi_stats

    from studio.app.optinist.core.edit_ROI.edit_ROI import CellType

    ops = ops.data
    stat = ops["stat"]
    temp_roi_data_dict = {**data.temp_add_roi, **data.temp_merge_roi}
    temp_roi_data_dict = dict(sorted(temp_roi_data_dict.items()))

    if "aspect" in ops:
        d0 = np.array([int(ops["aspect"] * 10), 10])
    else:
        d0 = ops["diameter"]
        if isinstance(d0, int):
            d0 = [d0, d0]

    for roi_id, roi_data in temp_roi_data_dict.items():
        iscell[int(roi_id)] = CellType.ROI
        # Save fluorescence traces of added roi
        if isinstance(roi_data, RoiPos):
            stat0 = get_stat0_add_roi(
                ops, roi_data.posx, roi_data.posy, roi_data.sizex, roi_data.sizey
            )
            Fcell, Fneu, _, _, _, ops, _ = masks_and_traces(ops, stat0, stat)
            ops["F"] = np.concatenate((ops["F"], Fcell), axis=0)
            ops["Fneu"] = np.concatenate((ops["Fneu"], Fneu), axis=0)

        # Save fluorescence traces of merged roi
        elif isinstance(roi_data, list):
            ypix = np.zeros((0,), np.int32)
            xpix = np.zeros((0,), np.int32)
            lam = np.zeros((0,), np.float32)

            merged_cells = np.unique(np.array(roi_data))
            F = np.mean((ops["F"][roi_data, :]), axis=0)
            Fneu = np.mean((ops["Fneu"][roi_data, :]), axis=0)

            # merge cells regions
            for n in merged_cells:
                ypix = np.append(ypix, stat[n]["ypix"])
                xpix = np.append(xpix, stat[n]["xpix"])
                lam = np.append(lam, stat[n]["lam"])

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

            ops["F"] = np.vstack((ops["F"], F))
            ops["Fneu"] = np.vstack((ops["Fneu"], Fneu))

            stat.append(stat0)
            stat = roi_stats(stat, d0[0], d0[1], ops["Ly"], ops["Lx"])
            stat[-1]["lam"] = stat[-1]["lam"] * merged_cells.size

    iscell[iscell == CellType.TEMP_DELETE] = CellType.NON_ROI

    ops["stat"] = stat
    data.commit()

    info = {
        "ops": Suite2pData(ops),
        "fluorescence": FluoData(ops["F"], file_name="fluorescence"),
        "iscell": IscellData(iscell),
        "cell_roi": RoiData(
            np.nanmax(data.im[iscell != CellType.NON_ROI], axis=0),
            output_dir=node_dirpath,
            file_name="cell_roi",
        ),
        "edit_roi_data": data,
        "nwbfile": set_nwbfile(ops, iscell, data, function_id),
    }
    return info
