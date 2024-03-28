from studio.app.common.dataclass import ImageData
from studio.app.optinist.core.nwb.nwb import NWBDATASET
from studio.app.optinist.dataclass import (
    EditRoiData,
    FluoData,
    IscellData,
    RoiData,
    Suite2pData,
)


def suite2p_roi(
    ops: Suite2pData, output_dir: str, params: dict = None, **kwargs
) -> dict(ops=Suite2pData, fluorescence=FluoData, iscell=IscellData):
    import numpy as np
    from suite2p import ROI, classification, default_ops, detection, extraction

    function_id = output_dir.split("/")[-1]
    print("start suite2p_roi:", function_id)

    nwbfile = kwargs.get("nwbfile", {})
    fs = nwbfile.get("imaging_plane", {}).get("imaging_rate", 30)

    ops = ops.data
    ops = {**default_ops(), **ops, **params, "fs": fs}

    # ROI detection
    ops_classfile = ops.get("classifier_path")
    builtin_classfile = classification.builtin_classfile
    user_classfile = classification.user_classfile
    if ops_classfile:
        print(f"NOTE: applying classifier {str(ops_classfile)}")
        classfile = ops_classfile
    elif ops["use_builtin_classifier"] or not user_classfile.is_file():
        print(f"NOTE: Applying builtin classifier at {str(builtin_classfile)}")
        classfile = builtin_classfile
    else:
        print(f"NOTE: applying default {str(user_classfile)}")
        classfile = user_classfile

    ops, stat = detection.detect(ops=ops, classfile=classfile)

    # ROI EXTRACTION
    ops, stat, F, Fneu, _, _ = extraction.create_masks_and_extract(ops, stat)
    stat = stat.tolist()

    # ROI CLASSIFICATION
    iscell = classification.classify(stat=stat, classfile=classfile)
    iscell = iscell[:, 0].astype(int)

    arrays = []
    for i, s in enumerate(stat):
        array = ROI(
            ypix=s["ypix"], xpix=s["xpix"], lam=s["lam"], med=s["med"], do_crop=False
        ).to_array(Ly=ops["Ly"], Lx=ops["Lx"])
        array *= i + 1
        arrays.append(array)

    im = np.stack(arrays)
    im[im == 0] = np.nan
    im -= 1

    # roiを追加
    roi_list = []
    for i in range(len(stat)):
        kargs = {}
        kargs["pixel_mask"] = np.array(
            [stat[i]["ypix"], stat[i]["xpix"], stat[i]["lam"]]
        ).T
        roi_list.append(kargs)

    # NWBを追加
    nwbfile = {}

    nwbfile[NWBDATASET.ROI] = {function_id: roi_list}
    nwbfile[NWBDATASET.POSTPROCESS] = {function_id: {"all_roi_img": im}}

    # iscellを追加
    nwbfile[NWBDATASET.COLUMN] = {
        function_id: {
            "name": "iscell",
            "description": "two columns - iscell & probcell",
            "data": iscell,
        }
    }

    # Fluorenceを追加
    nwbfile[NWBDATASET.FLUORESCENCE] = {
        function_id: {
            "Fluorescence": {
                "table_name": "Fluorescence",
                "region": list(range(len(F))),
                "name": "Fluorescence",
                "data": np.transpose(F),
                "unit": "lumens",
                "rate": ops["fs"],
            },
            "Neuropil": {
                "table_name": "Neuropil",
                "region": list(range(len(Fneu))),
                "name": "Neuropil",
                "data": np.transpose(Fneu),
                "unit": "lumens",
                "rate": ops["fs"],
            },
        }
    }

    ops["stat"] = stat
    ops["F"] = F
    ops["Fneu"] = Fneu

    info = {
        "ops": Suite2pData(ops),
        "max_proj": ImageData(
            ops["max_proj"], output_dir=output_dir, file_name="max_proj"
        ),
        "Vcorr": ImageData(ops["Vcorr"], output_dir=output_dir, file_name="Vcorr"),
        "fluorescence": FluoData(F, file_name="fluorescence"),
        "iscell": IscellData(iscell, file_name="iscell"),
        "all_roi": RoiData(
            np.nanmax(im, axis=0), output_dir=output_dir, file_name="all_roi"
        ),
        "non_cell_roi": RoiData(
            np.nanmax(im[iscell == 0], axis=0),
            output_dir=output_dir,
            file_name="noncell_roi",
        ),
        "cell_roi": RoiData(
            np.nanmax(im[iscell != 0], axis=0),
            output_dir=output_dir,
            file_name="cell_roi",
        ),
        "edit_roi_data": EditRoiData(images=ImageData(ops["filelist"]).data, im=im),
        "nwbfile": nwbfile,
    }

    return info
