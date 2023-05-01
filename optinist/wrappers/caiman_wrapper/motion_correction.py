from optinist.api.dataclass.dataclass import ImageData, RoiData
from optinist.api.nwb.nwb import NWBDATASET


def caiman_mc(
    image: ImageData, output_dir: str, params: dict = None
) -> dict(mc_images=ImageData):
    import numpy as np
    from caiman import load_memmap, save_memmap, stop_server
    from caiman.base.rois import extract_binary_masks_from_structural_channel
    from caiman.cluster import setup_cluster
    from caiman.motion_correction import MotionCorrect
    from caiman.source_extraction.cnmf.params import CNMFParams

    opts = CNMFParams()

    if params is not None:
        opts.change_params(params_dict=params)

    c, dview, n_processes = setup_cluster(
        backend="local", n_processes=None, single_thread=True
    )

    mc = MotionCorrect(image.path, dview=dview, **opts.get_group("motion"))

    mc.motion_correct(save_movie=True)
    border_to_0 = 0 if mc.border_nan == "copy" else mc.border_to_0

    # memory mapping
    fname_new = save_memmap(
        mc.mmap_file, base_name="memmap_", order="C", border_to_0=border_to_0
    )

    stop_server(dview=dview)

    # now load the file
    Yr, dims, T = load_memmap(fname_new)

    images = np.array(Yr.T.reshape((T,) + dims, order="F"))

    meanImg = images.mean(axis=0)
    rois = (
        extract_binary_masks_from_structural_channel(
            meanImg, gSig=7, expand_method="dilation"
        )[0]
        .reshape(meanImg.shape[0], meanImg.shape[1], -1)
        .transpose(2, 0, 1)
    )

    rois = rois.astype(np.float)

    for i, _ in enumerate(rois):
        rois[i] *= i + 1

    rois = np.nanmax(rois, axis=0)
    rois[rois == 0] = np.nan

    xy_trans_data = (
        (np.array(mc.x_shifts_els), np.array(mc.y_shifts_els))
        if params["pw_rigid"]
        else np.array(mc.shifts_rig)
    )

    mc_images = ImageData(images, output_dir=output_dir, file_name="mc_images")

    nwbfile = {}
    nwbfile[NWBDATASET.MOTION_CORRECTION] = {
        "caiman_mc": {
            "mc_data": mc_images,
            "xy_trans_data": xy_trans_data,
        }
    }

    info = {
        "mc_images": mc_images,
        "meanImg": ImageData(meanImg, output_dir=output_dir, file_name="meanImg"),
        "rois": RoiData(rois, output_dir=output_dir, file_name="rois"),
        "nwbfile": nwbfile,
    }

    return info
