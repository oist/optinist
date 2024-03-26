import gc

import numpy as np

from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.dataclass import ImageData
from studio.app.optinist.core.nwb.nwb import NWBDATASET
from studio.app.optinist.dataclass import EditRoiData, FluoData, IscellData, RoiData


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
        indx = np.argsort(patch_data)[::-1]

        if thr_method == "nrg":
            cumEn = np.cumsum(patch_data[indx] ** 2)
            if len(cumEn) == 0:
                pars = dict(
                    coordinates=np.array([]),
                    CoM=np.array([np.NaN, np.NaN]),
                    neuron_id=i + 1,
                )
                coordinates.append(pars)
                continue
            else:
                # we work with normalized values
                cumEn /= cumEn[-1]
                Bvec = np.ones(d)
                # we put it in a similar matrix
                Bvec[A.indices[A.indptr[i] : A.indptr[i + 1]][indx]] = cumEn
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


def util_recursive_flatten_params(params, result_params: dict, nest_counter=0):
    """
    Recursively flatten node parameters (operation for CaImAn CNMFParams)
    """
    # avoid infinite loops
    assert nest_counter <= 2, f"Nest depth overflow. [{nest_counter}]"
    nest_counter += 1

    for key, nested_param in params.items():
        if type(nested_param) is dict:
            util_recursive_flatten_params(nested_param, result_params, nest_counter)
        else:
            result_params[key] = nested_param


def caiman_cnmf(
    images: ImageData, output_dir: str, params: dict = None, **kwargs
) -> dict(fluorescence=FluoData, iscell=IscellData):
    from caiman import local_correlations, stop_server
    from caiman.cluster import setup_cluster
    from caiman.mmapping import prepare_shape
    from caiman.paths import memmap_frames_filename
    from caiman.source_extraction.cnmf import cnmf
    from caiman.source_extraction.cnmf.params import CNMFParams

    function_id = output_dir.split("/")[-1]
    print("start caiman_cnmf:", function_id)

    # flatten cmnf params segments.
    reshaped_params = {}
    util_recursive_flatten_params(params, reshaped_params)

    Ain = reshaped_params.pop("Ain", None)
    do_refit = reshaped_params.pop("do_refit", None)
    roi_thr = reshaped_params.pop("roi_thr", None)

    file_path = images.path
    if isinstance(file_path, list):
        file_path = file_path[0]

    images = images.data

    # np.arrayをmmapへ変換
    order = "C"
    dims = images.shape[1:]
    T = images.shape[0]
    shape_mov = (np.prod(dims), T)

    dir_path = join_filepath(file_path.split("/")[:-1])
    basename = file_path.split("/")[-1]
    fname_tot = memmap_frames_filename(basename, dims, T, order)

    mmap_images = np.memmap(
        join_filepath([dir_path, fname_tot]),
        mode="w+",
        dtype=np.float32,
        shape=prepare_shape(shape_mov),
        order=order,
    )

    mmap_images = np.reshape(mmap_images.T, [T] + list(dims), order="F")
    mmap_images[:] = images[:]

    del images
    gc.collect()

    nwbfile = kwargs.get("nwbfile", {})
    fr = nwbfile.get("imaging_plane", {}).get("imaging_rate", 30)

    if reshaped_params is None:
        ops = CNMFParams()
    else:
        ops = CNMFParams(params_dict={**reshaped_params, "fr": fr})

    if "dview" in locals():
        stop_server(dview=dview)  # noqa: F821

    c, dview, n_processes = setup_cluster(
        backend="local", n_processes=None, single_thread=True
    )

    cnm = cnmf.CNMF(n_processes=n_processes, dview=dview, Ain=Ain, params=ops)
    cnm = cnm.fit(mmap_images)

    if do_refit:
        cnm = cnm.refit(mmap_images, dview=dview)

    cnm.estimates.evaluate_components(mmap_images, cnm.params, dview=dview)

    stop_server(dview=dview)

    # contours plot
    Cn = local_correlations(mmap_images.transpose(1, 2, 0))
    Cn[np.isnan(Cn)] = 0

    thr_method = "nrg"
    swap_dim = False

    idx_good = cnm.estimates.idx_components
    idx_bad = cnm.estimates.idx_components_bad
    if not isinstance(idx_good, list):
        idx_good = idx_good.tolist()
    if not isinstance(idx_bad, list):
        idx_bad = idx_bad.tolist()

    iscell = np.concatenate([np.ones(len(idx_good)), np.zeros(len(idx_bad))])

    cell_ims = get_roi(
        cnm.estimates.A[:, idx_good], roi_thr, thr_method, swap_dim, dims
    )
    cell_ims = np.stack(cell_ims).astype(float)
    cell_ims[cell_ims == 0] = np.nan
    cell_ims -= 1
    n_rois = len(cell_ims)

    if len(idx_bad) > 0:
        non_cell_ims = get_roi(
            cnm.estimates.A[:, idx_bad], roi_thr, thr_method, swap_dim, dims
        )
        non_cell_ims = np.stack(non_cell_ims).astype(float)
        for i, j in enumerate(range(n_rois, n_rois + len(non_cell_ims))):
            non_cell_ims[i, :] = np.where(non_cell_ims[i, :] != 0, j, 0)
        non_cell_roi = np.nanmax(non_cell_ims, axis=0).astype(float)
    else:
        non_cell_ims = np.zeros((0, *dims))
        non_cell_roi = np.zeros(dims)
        non_cell_roi[non_cell_roi == 0] = np.nan
    non_cell_ims[non_cell_ims == 0] = np.nan

    n_noncell_rois = len(non_cell_ims)

    im = np.vstack([cell_ims, non_cell_ims])

    # NWBの追加
    nwbfile = {}
    # NWBにROIを追加
    roi_list = []
    n_cells = cnm.estimates.A.shape[-1]
    for i in range(n_cells):
        kargs = {}
        kargs["image_mask"] = cnm.estimates.A.T[i].T.toarray().reshape(dims)
        if hasattr(cnm.estimates, "accepted_list"):
            kargs["accepted"] = i in cnm.estimates.accepted_list
        if hasattr(cnm.estimates, "rejected_list"):
            kargs["rejected"] = i in cnm.estimates.rejected_list
        roi_list.append(kargs)

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

    # Fluorescence
    fluorescence = cnm.estimates.C

    nwbfile[NWBDATASET.FLUORESCENCE] = {
        function_id: {
            "Fluorescence": {
                "table_name": "ROIs",
                "region": list(range(n_rois + n_noncell_rois)),
                "name": "Fluorescence",
                "data": fluorescence.T,
                "unit": "lumens",
            }
        }
    }

    info = {
        "images": ImageData(
            np.array(Cn * 255, dtype=np.uint8),
            output_dir=output_dir,
            file_name="images",
        ),
        "fluorescence": FluoData(fluorescence, file_name="fluorescence"),
        "iscell": IscellData(iscell, file_name="iscell"),
        "all_roi": RoiData(
            np.nanmax(im, axis=0), output_dir=output_dir, file_name="all_roi"
        ),
        "cell_roi": RoiData(
            np.nanmax(im[iscell != 0], axis=0),
            output_dir=output_dir,
            file_name="cell_roi",
        ),
        "non_cell_roi": RoiData(
            non_cell_roi, output_dir=output_dir, file_name="non_cell_roi"
        ),
        "edit_roi_data": EditRoiData(mmap_images, im),
        "nwbfile": nwbfile,
    }

    return info
