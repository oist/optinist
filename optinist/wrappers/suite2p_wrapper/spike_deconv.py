from optinist.api.dataclass.dataclass import FluoData, Suite2pData
from optinist.api.nwb.nwb import NWBDATASET


def suite2p_spike_deconv(
    ops: Suite2pData, output_dir: str, params: dict = None
) -> dict(ops=Suite2pData, spks=FluoData):
    import numpy as np
    from suite2p import default_ops, extraction

    print("start suite2_spike_deconv")
    ops = ops.data

    ops = {**default_ops(), **ops, **params}

    dF = ops["F"].copy() - ops["neucoeff"] * ops["Fneu"]
    dF = extraction.preprocess(
        F=dF,
        baseline=ops["baseline"],
        win_baseline=ops["win_baseline"],
        sig_baseline=ops["sig_baseline"],
        fs=ops["fs"],
        prctile_baseline=ops["prctile_baseline"],
    )
    spks = extraction.oasis(
        F=dF, batch_size=ops["batch_size"], tau=ops["tau"], fs=ops["fs"]
    )

    ops["spks"] = spks

    # NWBを追加
    nwbfile = {}

    # roiを追加
    stat = ops["stat"]
    roi_list = []
    for i in range(len(stat)):
        kargs = {}
        kargs["pixel_mask"] = np.array(
            [stat[i]["ypix"], stat[i]["xpix"], stat[i]["lam"]]
        ).T
        roi_list.append(kargs)

    nwbfile[NWBDATASET.ROI] = {"roi_list": roi_list}

    # Fluorenceを追加
    name = "Deconvolved"
    nwbfile[NWBDATASET.FLUORESCENCE] = {
        name: {
            "table_name": name,
            "region": list(range(len(spks))),
            "name": name,
            "data": spks,
            "unit": "lumens",
            "rate": ops["fs"],
        }
    }

    info = {
        "ops": Suite2pData(ops),
        "spks": FluoData(spks, file_name="spks"),
        "nwbfile": nwbfile,
    }

    return info
