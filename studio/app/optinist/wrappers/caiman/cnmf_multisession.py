import gc

import imageio
import numpy as np

from studio.app.common.core.experiment.experiment import ExptOutputPathIds
from studio.app.common.core.logger import AppLogger
from studio.app.common.dataclass import ImageData
from studio.app.optinist.core.nwb.nwb import NWBDATASET
from studio.app.optinist.dataclass import EditRoiData, FluoData, IscellData, RoiData
from studio.app.optinist.wrappers.caiman.cnmf import (
    get_roi,
    nwb_image_mask,
    util_cleanup_image_memmap,
    util_download_model_files,
    util_get_image_memmap,
)
from studio.app.optinist.wrappers.optinist.utils import recursive_flatten_params

logger = AppLogger.get_logger()


def caiman_cnmf_multisession(
    images: ImageData, output_dir: str, params: dict = None, **kwargs
) -> dict(fluorescence=FluoData, iscell=IscellData):
    from caiman import load, local_correlations, stop_server
    from caiman.base.rois import register_multisession
    from caiman.cluster import setup_cluster
    from caiman.source_extraction.cnmf import cnmf
    from caiman.source_extraction.cnmf.params import CNMFParams

    function_id = ExptOutputPathIds(output_dir).function_id
    logger.info(f"start caiman_cnmf_multisession: {function_id}")

    # NOTE: evaluate_components requires cnn_model files in caiman_data directory.
    util_download_model_files()

    flattened_params = {}
    recursive_flatten_params(params, flattened_params)
    params = flattened_params

    Ain = params.pop("Ain", None)
    roi_thr = params.pop("roi_thr", None)

    # mulisiession params
    n_reg_files = params.pop("n_reg_files", 2)
    if n_reg_files < 2:
        raise Exception(f"Set n_reg_files to a integer value gte 2. Now {n_reg_files}.")
    reg_file_rate = params.pop("reg_file_rate", 1.0)
    if reg_file_rate > 1.0:
        logger.warn(
            f"reg_file_rate {reg_file_rate}, should be lte 1. Using 1.0 instead."
        )
        reg_file_rate = 1.0

    align_flag = params.pop("align_flag", True)
    max_thr = params.pop("max_thr", 0)
    use_opt_flow = params.pop("use_opt_flow", True)
    thresh_cost = params.pop("thresh_cost", 0.7)
    max_dist = params.pop("max_dist", 10)
    enclosed_thr = params.pop("enclosed_thr", None)

    split_image_paths = images.split_image(output_dir, n_files=n_reg_files)
    n_split_images = len(split_image_paths)

    logger.info(f"image was split into {n_split_images} parts.")

    nwbfile = kwargs.get("nwbfile", {})
    fr = nwbfile.get("imaging_plane", {}).get("imaging_rate", 30)

    if params is None:
        ops = CNMFParams()
    else:
        ops = CNMFParams(params_dict={**params, "fr": fr})

    if "dview" in locals():
        stop_server(dview=dview)  # noqa: F821

    # TODO: Add parameters for node
    n_processes = 1
    dview = None
    # This process launches another process to run the CNMF algorithm,
    # so this node use at least 2 core.
    if n_processes == 1:
        c, dview, n_processes = setup_cluster(
            backend="single", n_processes=n_processes, single_thread=True
        )
    else:
        c, dview, n_processes = setup_cluster(
            backend="multiprocessing", n_processes=n_processes
        )
    logger.info(f"n_processes: {n_processes}")

    cnm_list = []
    templates = []
    mmap_paths = []
    for split_image_path in split_image_paths:
        split_image = imageio.volread(split_image_path)
        split_image_mmap, _, mmap_path = util_get_image_memmap(
            function_id, split_image, split_image_path
        )
        mmap_paths.append(mmap_path)
        del split_image
        gc.collect()

        # ops.change_params("fnames", [image_path])
        cnm = cnmf.CNMF(n_processes=n_processes, dview=dview, Ain=Ain, params=ops)
        cnm = cnm.fit(split_image_mmap)
        cnm_list.append(cnm)
        templates.append(load(split_image_path).mean(0))

        del split_image_mmap
        gc.collect()

    stop_server(dview=dview)

    spatial = [cnm.estimates.A for cnm in cnm_list]
    dims = templates[0].shape

    spatial_union, assignments, matchings = register_multisession(
        A=spatial,
        dims=dims,
        templates=templates,
        align_flag=align_flag,
        max_thr=max_thr,
        use_opt_flow=use_opt_flow,
        thresh_cost=thresh_cost,
        max_dist=max_dist,
        enclosed_thr=enclosed_thr,
    )

    reg_files = int(len(split_image_paths) * reg_file_rate)
    assignments_filtered = np.array(
        np.nan_to_num(
            assignments[np.sum(~np.isnan(assignments), axis=1) >= reg_files],
        ),
        dtype=int,
    )
    if assignments_filtered.shape[0] == 0:
        raise Exception(
            f"No same cells found in {reg_files} or more "
            + f"out of {n_split_images} files. "
            + "Try to set lower reg_file_rate.",
        )

    spatial_filtered = spatial[0][:, assignments_filtered[:, 0]]

    fluorescence = np.array(
        [
            np.concatenate(
                [
                    cnm_list[j].estimates.C[int(assignments_filtered[i, j])]
                    for j in range(assignments_filtered.shape[1])
                ]
            )
            for i in range(assignments_filtered.shape[0])
        ]
    )

    # contours plot
    thr_method = "nrg"
    swap_dim = False

    iscell = np.concatenate([np.ones(assignments_filtered.shape[0], dtype=int)])

    cell_ims = get_roi(spatial_filtered, roi_thr, thr_method, swap_dim, dims)
    cell_ims = np.stack(cell_ims).astype(float)
    cell_ims[cell_ims == 0] = np.nan
    cell_ims -= 1
    n_rois = len(cell_ims)

    # NWBの追加
    nwbfile = {}
    # NWBにROIを追加
    roi_list = []
    n_cells = spatial_filtered.shape[-1]
    for i in range(n_cells):
        kargs = {}
        kargs["image_mask"] = nwb_image_mask(spatial_filtered, i, dims, swap_dim)
        roi_list.append(kargs)

    nwbfile[NWBDATASET.ROI] = {function_id: {"roi_list": roi_list}}
    nwbfile[NWBDATASET.POSTPROCESS] = {function_id: {"all_roi_img": cell_ims}}

    # iscellを追加
    nwbfile[NWBDATASET.COLUMN] = {
        function_id: {
            "name": "iscell",
            "description": "iscell",
            "data": iscell,
        }
    }

    nwbfile[NWBDATASET.FLUORESCENCE] = {
        function_id: {
            "Fluorescence": {
                "table_name": "ROIs",
                "region": list(range(n_rois)),
                "name": "Fluorescence",
                "data": fluorescence,
                "unit": "lumens",
                "rate": fr,
                "comments": f"ROIs were detected in {reg_files} out of "
                + f"{n_split_images} sessions ({reg_file_rate} overall).",
            }
        }
    }

    # get mean image
    file_path = images.path
    if isinstance(file_path, list):
        file_path = file_path[0]
    images = images.data
    mmap_images, dims, mmap_path = util_get_image_memmap(
        function_id, images.data, file_path
    )
    mmap_paths.append(mmap_path)

    Cn = local_correlations(mmap_images.transpose(1, 2, 0))
    Cn[np.isnan(Cn)] = 0

    info = {
        "images": ImageData(
            np.array(Cn * 255, dtype=np.uint8),
            output_dir=output_dir,
            file_name="images",
        ),
        "fluorescence": FluoData(fluorescence, file_name="fluorescence"),
        "iscell": IscellData(iscell, file_name="iscell"),
        "cell_roi": RoiData(
            np.nanmax(cell_ims[iscell != 0], axis=0),
            output_dir=output_dir,
            file_name="cell_roi",
        ),
        "edit_roi_data": EditRoiData(mmap_images, cell_ims),
        "nwbfile": nwbfile,
    }

    # Clean up temporary files
    try:
        util_cleanup_image_memmap(mmap_paths)
    except Exception as e:
        logger.error("Failed to cleanup memmap files.")
        logger.error(e)

    return info
