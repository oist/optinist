import gc
import os

import numpy as np
import requests

from studio.app.common.core.experiment.experiment import ExptOutputPathIds
from studio.app.common.core.logger import AppLogger
from studio.app.common.core.utils.filepath_creater import (
    create_directory,
    join_filepath,
)
from studio.app.common.dataclass import ImageData
from studio.app.optinist.core.nwb.nwb import NWBDATASET
from studio.app.optinist.dataclass import EditRoiData, FluoData, IscellData, RoiData
from studio.app.optinist.wrappers.optinist.utils import recursive_flatten_params

logger = AppLogger.get_logger()


def get_roi(A, roi_thr, thr_method, swap_dim, dims):
    from scipy.ndimage import binary_fill_holes
    from skimage.measure import find_contours

    d, nr = np.shape(A)

    # for each patches
    ims = []
    coordinates = []
    for i in range(nr):
        pars = dict()
        # we compute the cumulative sum of the energy of the Ath component
        # that has been ordered from least to highest
        patch_data = A.data[A.indptr[i] : A.indptr[i + 1]]
        idx = np.argsort(patch_data)[::-1]

        if thr_method == "nrg":
            cumEn = np.cumsum(patch_data[idx] ** 2)
            if len(cumEn) == 0:
                pars = dict(
                    coordinates=np.array([]),
                    CoM=np.array([np.NaN, np.NaN]),
                    neuron_id=i + 1,
                )
                coordinates.append(pars)
                # keep one mask per column so row i stays component i
                ims.append(np.zeros(dims))
                continue
            else:
                # we work with normalized values
                cumEn /= cumEn[-1]
                Bvec = np.ones(d)
                # we put it in a similar matrix
                Bvec[A.indices[A.indptr[i] : A.indptr[i + 1]][idx]] = cumEn
        else:
            Bvec = np.zeros(d)
            Bvec[A.indices[A.indptr[i] : A.indptr[i + 1]]] = (
                patch_data / patch_data.max()
            )

        if swap_dim:
            Bmat = np.reshape(Bvec, dims, order="C")
        else:
            Bmat = np.reshape(Bvec, dims, order="F")

        r_mask = np.zeros_like(Bmat, dtype="bool")
        contour = find_contours(Bmat, roi_thr)
        for c in contour:
            r_mask[np.round(c[:, 0]).astype("int"), np.round(c[:, 1]).astype("int")] = 1

        # Fill in the hole created by the contour boundary
        r_mask = binary_fill_holes(r_mask)
        ims.append(r_mask + (i * r_mask))

    return ims


def util_get_image_memmap(function_id: str, images: np.ndarray, file_path: str):
    """
    convert np.ndarray to mmap
    """
    from caiman.mmapping import prepare_shape
    from caiman.paths import memmap_frames_filename

    order = "C"
    dims = images.shape[1:]
    T = images.shape[0]
    shape_mov = (np.prod(dims), T)

    dir_path = join_filepath(file_path.split("/")[:-1])
    file_basename = file_path.split("/")[-1]
    mmap_basename = f"{file_basename}.{function_id}"
    fname_tot = memmap_frames_filename(mmap_basename, dims, T, order)
    mmap_path = join_filepath([dir_path, fname_tot])

    mmap_images = np.memmap(
        mmap_path,
        mode="w+",
        dtype=np.float32,
        shape=prepare_shape(shape_mov),
        order=order,
    )

    mmap_images = np.reshape(mmap_images.T, [T] + list(dims), order="F")
    mmap_images[:] = images[:]
    return mmap_images, dims, mmap_path


def util_cleanup_image_memmap(mmap_paths: list):
    for mmap_path in mmap_paths:
        if mmap_path.endswith(".mmap") and os.path.isfile(mmap_path):
            os.remove(mmap_path)


def util_download_model_files():
    """
    download model files for component evaluation
    """
    # NOTE: We specify the version of the CaImAn to download.
    base_url = "https://github.com/flatironinstitute/CaImAn/raw/v1.9.12/model"
    model_files = [
        "cnn_model.h5",
        "cnn_model.h5.pb",
        "cnn_model.json",
        "cnn_model_online.h5",
        "cnn_model_online.h5.pb",
        "cnn_model_online.json",
    ]

    caiman_data_dir = os.path.join(os.path.expanduser("~"), "caiman_data")
    if not os.path.exists(caiman_data_dir):
        create_directory(caiman_data_dir)

    model_dir = join_filepath([caiman_data_dir, "model"])
    if not os.path.exists(model_dir):
        create_directory(join_filepath(model_dir))

    if len(os.listdir(model_dir)) < len(model_files):
        for model in model_files:
            url = f"{base_url}/{model}"
            file_path = join_filepath([model_dir, model])
            if not os.path.exists(file_path):
                logger.info(f"Downloading {model}")
                response = requests.get(url)
                with open(file_path, "wb") as f:
                    f.write(response.content)


def nwb_image_mask(A, comp, dims, swap_dim):
    # the pixel order get_roi reads A in, so the NWB masks match the ROI images
    return A[:, comp].toarray().reshape(dims, order="C" if swap_dim else "F")


def component_outputs(
    A,
    C,
    idx_good,
    idx_bad,
    dims,
    roi_thr,
    thr_method,
    swap_dim,
    accepted_list=None,
    rejected_list=None,
    n_frames=0,
):
    """im, iscell, the NWB ROI table and fluorescence in one shared row order.

    Row i of every output is component order[i], accepted components first, so
    the ROI editor and the analyses can index fluorescence by im row.
    """
    order = list(idx_good) + list(idx_bad)
    n_components = len(order)

    iscell = np.concatenate(
        [
            np.ones(len(idx_good), dtype=int),
            np.zeros(len(idx_bad), dtype=int),
        ]
    )

    if len(idx_good) > 0 and A is not None:
        cell_ims = get_roi(A[:, idx_good], roi_thr, thr_method, swap_dim, dims)
        if len(cell_ims) != len(idx_good):
            raise ValueError("one mask per accepted component expected")
        if len(cell_ims) > 0:  # Check if get_roi returned any ROIs
            cell_ims = np.stack(cell_ims).astype(float)
            cell_ims[cell_ims == 0] = np.nan
            cell_ims = np.where(
                np.isnan(cell_ims), cell_ims, cell_ims - 1
            )  # Safer subtraction
            n_rois = len(cell_ims)
        else:
            cell_ims = np.zeros((0, *dims))
            n_rois = 0
    else:
        cell_ims = np.zeros((0, *dims))
        n_rois = 0

    if len(idx_bad) > 0:
        non_cell_ims = get_roi(A[:, idx_bad], roi_thr, thr_method, swap_dim, dims)
        if len(non_cell_ims) != len(idx_bad):
            raise ValueError("one mask per rejected component expected")
        non_cell_ims = np.stack(non_cell_ims).astype(float)
        for i, j in enumerate(range(n_rois, n_rois + len(non_cell_ims))):
            non_cell_ims[i, :] = np.where(non_cell_ims[i, :] != 0, j + 1, 0)
    else:
        non_cell_ims = np.zeros((0, *dims))
    non_cell_ims[non_cell_ims == 0] = np.nan
    non_cell_ims -= 1

    n_noncell_rois = len(non_cell_ims)
    non_cell_roi = (
        np.nanmax(non_cell_ims, axis=0) if n_noncell_rois > 0 else np.full(dims, np.nan)
    )

    im = (
        np.vstack([cell_ims, non_cell_ims])
        if n_components > 0
        else np.zeros((0, *dims))
    )

    roi_list = []
    for comp in order:
        kargs = {}
        kargs["image_mask"] = nwb_image_mask(A, comp, dims, swap_dim)
        if accepted_list is not None:
            kargs["accepted"] = comp in accepted_list
        if rejected_list is not None:
            kargs["rejected"] = comp in rejected_list
        roi_list.append(kargs)

    fluorescence = (
        np.asarray(C)[order]
        if n_components > 0 and C is not None
        else np.zeros((0, n_frames))
    )

    return im, iscell, roi_list, fluorescence, non_cell_roi, n_rois, n_noncell_rois


def caiman_cnmf(
    images: ImageData, output_dir: str, params: dict = None, **kwargs
) -> dict(fluorescence=FluoData, iscell=IscellData):
    from caiman import local_correlations, stop_server
    from caiman.cluster import setup_cluster
    from caiman.source_extraction.cnmf import cnmf, online_cnmf
    from caiman.source_extraction.cnmf.params import CNMFParams

    function_id = ExptOutputPathIds(output_dir).function_id
    logger.info(f"start caiman_cnmf: {function_id}")

    # NOTE: evaluate_components requires cnn_model files in caiman_data directory.
    util_download_model_files()

    flattened_params = {}
    recursive_flatten_params(params, flattened_params)
    params = flattened_params

    Ain = params.pop("Ain", None)
    do_refit = params.pop("do_refit", None)
    roi_thr = params.pop("roi_thr", None)
    use_online = params.pop("use_online", False)

    file_path = images.path
    if isinstance(file_path, list):
        file_path = file_path[0]

    images = images.data
    mmap_paths = []
    mmap_images, dims, mmap_path = util_get_image_memmap(function_id, images, file_path)
    mmap_paths.append(mmap_path)

    del images
    gc.collect()

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

    if use_online:
        ops.change_params(
            {
                "fnames": [mmap_path],
                # NOTE: These params uses np.inf as default in CaImAn.
                # Yaml cannot serialize np.inf, so default value in yaml is None.
                "max_comp_update_shape": params["max_comp_update_shape"] or np.inf,
                "num_times_comp_updated": params["update_num_comps"] or np.inf,
            }
        )
        cnm = online_cnmf.OnACID(dview=dview, Ain=Ain, params=ops)
        cnm.fit_online()
    else:
        cnm = cnmf.CNMF(n_processes=n_processes, dview=dview, Ain=Ain, params=ops)
        cnm = cnm.fit(mmap_images)

        if do_refit:
            cnm = cnm.refit(mmap_images, dview=dview)

    # Check if any components were found
    n_components = cnm.estimates.A.shape[1] if hasattr(cnm.estimates, "A") else 0

    if n_components > 0:
        # Only evaluate components if we found some
        cnm.estimates.evaluate_components(mmap_images, cnm.params, dview=dview)
        idx_good = cnm.estimates.idx_components
        idx_bad = cnm.estimates.idx_components_bad
        if not isinstance(idx_good, list):
            idx_good = idx_good.tolist()
        if not isinstance(idx_bad, list):
            idx_bad = idx_bad.tolist()
    else:
        # No components found
        idx_good = []
        idx_bad = []

    stop_server(dview=dview)

    # contours plot
    Cn = local_correlations(mmap_images.transpose(1, 2, 0))
    Cn[np.isnan(Cn)] = 0

    thr_method = "nrg"
    swap_dim = False

    (
        im,
        iscell,
        roi_list,
        fluorescence,
        non_cell_roi,
        n_rois,
        n_noncell_rois,
    ) = component_outputs(
        getattr(cnm.estimates, "A", None),
        getattr(cnm.estimates, "C", None),
        idx_good,
        idx_bad,
        dims,
        roi_thr,
        thr_method,
        swap_dim,
        getattr(cnm.estimates, "accepted_list", None),
        getattr(cnm.estimates, "rejected_list", None),
        mmap_images.shape[0],
    )

    # NWB additions
    nwbfile = {}
    nwbfile[NWBDATASET.ROI] = {function_id: {"roi_list": roi_list}}
    nwbfile[NWBDATASET.POSTPROCESS] = {function_id: {"all_roi_img": im}}

    # Add iscell to NWB
    nwbfile[NWBDATASET.COLUMN] = {
        function_id: {
            "name": "iscell",
            "description": "two columns - iscell & probcell",
            "data": iscell,
        }
    }

    nwbfile[NWBDATASET.FLUORESCENCE] = {
        function_id: {
            "Fluorescence": {
                "table_name": "ROIs",
                "region": list(range(n_rois + n_noncell_rois)),
                "name": "Fluorescence",
                "data": fluorescence,
                "unit": "lumens",
                "rate": fr,
            }
        }
    }

    # Create empty ROI array filled with NaN for when no components found
    empty_roi = np.full(dims, np.nan)

    info = {
        "images": ImageData(
            np.array(Cn * 255, dtype=np.uint8),
            output_dir=output_dir,
            file_name="images",
        ),
        "fluorescence": FluoData(fluorescence, file_name="fluorescence"),
        "iscell": IscellData(iscell, file_name="iscell"),
        "all_roi": RoiData(
            np.nanmax(im, axis=0) if len(im) > 0 else empty_roi.copy(),
            output_dir=output_dir,
            file_name="all_roi",
        ),
        "cell_roi": RoiData(
            np.nanmax(im[iscell != 0], axis=0, initial=np.nan)
            if len(im) > 0
            else empty_roi.copy(),
            output_dir=output_dir,
            file_name="cell_roi",
        ),
        "non_cell_roi": RoiData(
            non_cell_roi if len(im) > 0 else empty_roi.copy(),
            output_dir=output_dir,
            file_name="non_cell_roi",
        ),
        "edit_roi_data": EditRoiData(mmap_images, im),
        "nwbfile": nwbfile,
    }

    # Clean up temporary files
    try:
        util_cleanup_image_memmap(mmap_paths)
    except Exception as e:
        logger.error("Failed to cleanup memmap files.")
        logger.error(e)

    return info
