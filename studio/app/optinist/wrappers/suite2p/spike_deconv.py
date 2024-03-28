from studio.app.optinist.core.nwb.nwb import NWBDATASET
from studio.app.optinist.dataclass import SpikingActivityData, Suite2pData


def suite2p_spike_deconv(
    ops: Suite2pData, output_dir: str, params: dict = None, **kwargs
) -> dict(ops=Suite2pData, spks=SpikingActivityData):
    import numpy as np
    from suite2p import default_ops, extraction

    function_id = output_dir.split("/")[-1]
    print("start suite2_spike_deconv:", function_id)

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

    nwbfile[NWBDATASET.ROI] = {function_id: roi_list}

    # Fluorenceを追加
    nwbfile[NWBDATASET.FLUORESCENCE] = {
        function_id: {
            "Deconvolved": {
                "table_name": "Deconvolved",
                "region": list(range(len(spks))),
                "name": function_id + "_Deconvolved",
                "data": spks,
                "unit": "lumens",
                "rate": ops["fs"],
            }
        }
    }

    info = {
        "ops": Suite2pData(ops),
        "spks": SpikingActivityData(spks, file_name="spks"),
        "nwbfile": nwbfile,
    }

    return info
