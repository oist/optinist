import gc

import numpy as np

from optinist.api.dataclass.dataclass import (
    CaimanCnmfData,
    FluoData,
    ImageData,
    IscellData,
    RoiData,
)
from optinist.api.nwb.nwb import NWBDATASET
from optinist.api.utils.filepath_creater import join_filepath


def get_roi(A, thr, thr_method, swap_dim, dims):
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
        contour = find_contours(Bmat, thr)
        for c in contour:
            r_mask[np.round(c[:, 0]).astype("int"), np.round(c[:, 1]).astype("int")] = 1

        # Fill in the hole created by the contour boundary
        r_mask = binary_fill_holes(r_mask)
        ims.append(r_mask + (i * r_mask))

    return ims


def caiman_cnmf(
    images: ImageData, output_dir: str, params: dict = None
) -> dict(fluorescence=FluoData, iscell=IscellData):
    import scipy
    from caiman import local_correlations, stop_server
    from caiman.cluster import setup_cluster
    from caiman.mmapping import prepare_shape
    from caiman.paths import memmap_frames_filename
    from caiman.source_extraction.cnmf import cnmf
    from caiman.source_extraction.cnmf.params import CNMFParams

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

    if params is None:
        ops = CNMFParams()
    else:
        ops = CNMFParams(params_dict=params)

    if "dview" in locals():
        stop_server(dview=dview)  # noqa: F821

    c, dview, n_processes = setup_cluster(
        backend="local", n_processes=None, single_thread=True
    )

    cnm = cnmf.CNMF(n_processes=n_processes, dview=dview, Ain=None, params=ops)
    cnm = cnm.fit(mmap_images)

    stop_server(dview=dview)

    # contours plot
    Cn = local_correlations(mmap_images.transpose(1, 2, 0))
    Cn[np.isnan(Cn)] = 0

    thr = params["thr"]
    thr_method = "nrg"
    swap_dim = False

    iscell = np.concatenate(
        [np.ones(cnm.estimates.A.shape[-1]), np.zeros(cnm.estimates.b.shape[-1])]
    )

    ims = get_roi(cnm.estimates.A, thr, thr_method, swap_dim, dims)
    ims = np.stack(ims)
    cell_roi = np.nanmax(ims, axis=0).astype(float)
    cell_roi[cell_roi == 0] = np.nan

    non_cell_roi_ims = get_roi(
        scipy.sparse.csc_matrix(cnm.estimates.b), thr, thr_method, swap_dim, dims
    )
    non_cell_roi_ims = np.stack(non_cell_roi_ims)
    non_cell_roi = np.nanmax(non_cell_roi_ims, axis=0).astype(float)
    non_cell_roi[non_cell_roi == 0] = np.nan

    all_roi = np.nanmax(np.stack([cell_roi, non_cell_roi]), axis=0)

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

    # backgroundsを追加
    bg_list = []
    for bg in cnm.estimates.b.T:
        kargs = {}
        kargs["image_mask"] = bg.reshape(dims)
        if hasattr(cnm.estimates, "accepted_list"):
            kargs["accepted"] = False
        if hasattr(cnm.estimates, "rejected_list"):
            kargs["rejected"] = False
        bg_list.append(kargs)

    nwbfile[NWBDATASET.ROI] = {
        "roi_list": roi_list,
        "bg_list": bg_list,
    }

    # iscellを追加
    nwbfile[NWBDATASET.COLUMN] = {
        "roi_column": {
            "name": "iscell",
            "discription": "two columns - iscell & probcell",
            "data": iscell,
        }
    }

    # Fluorescence
    n_rois = cnm.estimates.A.shape[-1]
    n_bg = len(cnm.estimates.f)

    nwbfile[NWBDATASET.FLUORESCENCE] = {
        "RoiResponseSeries": {
            "table_name": "ROIs",
            "region": list(range(n_rois)),
            "name": "RoiResponseSeries",
            "data": cnm.estimates.C.T,
            "unit": "lumens",
        },
        "Background_Fluorescence_Response": {
            "table_name": "Background",
            "region": list(range(n_rois, n_rois + n_bg)),
            "name": "Background_Fluorescence_Response",
            "data": cnm.estimates.f.T,
            "unit": "lumens",
        },
    }

    fluorescence = np.concatenate(
        [
            cnm.estimates.C,
            cnm.estimates.f,
        ]
    )

    cnmf_data = {}
    cnmf_data["fluorescence"] = fluorescence
    cnmf_data["im"] = np.concatenate([ims, non_cell_roi_ims], axis=0)
    cnmf_data["is_cell"] = iscell.astype(bool)
    cnmf_data["images"] = mmap_images

    info = {
        "images": ImageData(
            np.array(Cn * 255, dtype=np.uint8),
            output_dir=output_dir,
            file_name="images",
        ),
        "fluorescence": FluoData(fluorescence, file_name="fluorescence"),
        "iscell": IscellData(iscell, file_name="iscell"),
        "all_roi": RoiData(all_roi, output_dir=output_dir, file_name="all_roi"),
        "cell_roi": RoiData(cell_roi, output_dir=output_dir, file_name="cell_roi"),
        "non_cell_roi": RoiData(
            non_cell_roi, output_dir=output_dir, file_name="non_cell_roi"
        ),
        "nwbfile": nwbfile,
        "cnmf_data": CaimanCnmfData(cnmf_data),
    }

    return info
