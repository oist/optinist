from wrappers.data_wrapper import *
from cui_api.utils import join_file_path
import gc
from wrappers.nwb_wrapper.const import NWBDATASET


def get_roi(A, thr, thr_method, swap_dim, dims):
    from skimage.measure import find_contours
    from scipy.ndimage import binary_fill_holes

    d, nr = np.shape(A)

    # for each patches
    ims = []
    for i in range(nr):
        pars = dict()
        # we compute the cumulative sum of the energy of the Ath component that has been ordered from least to highest
        patch_data = A.data[A.indptr[i]:A.indptr[i + 1]]
        indx = np.argsort(patch_data)[::-1]

        if thr_method == 'nrg':
            cumEn = np.cumsum(patch_data[indx]**2)
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
                Bvec[A.indices[A.indptr[i]:A.indptr[i + 1]][indx]] = cumEn
        else:
            Bvec = np.zeros(d)
            Bvec[A.indices[A.indptr[i]:A.indptr[i + 1]]] = patch_data / patch_data.max()

        if swap_dim:
            Bmat = np.reshape(Bvec, dims, order='C')
        else:
            Bmat = np.reshape(Bvec, dims, order='F')

        r_mask = np.zeros_like(Bmat, dtype='bool')
        contour = find_contours(Bmat, thr)
        for c in contour:
            r_mask[np.round(c[:, 0]).astype('int'), np.round(c[:, 1]).astype('int')] = 1
        
        # Fill in the hole created by the contour boundary
        r_mask = binary_fill_holes(r_mask)
        ims.append(r_mask + (i * r_mask))

    return ims


def caiman_cnmf(
        images: ImageData,
        nwbfile: NWBFile=None,
        params: dict=None
    ) -> {'fluorescence': FluoData, 'iscell': IscellData}:
    import caiman
    from caiman import local_correlations, stop_server
    from caiman.paths import memmap_frames_filename
    from caiman.mmapping import prepare_shape
    from caiman.cluster import setup_cluster
    from caiman.source_extraction.cnmf import cnmf
    from caiman.source_extraction.cnmf.params import CNMFParams
    import caiman.utils.visualization as visualization
    import numpy as np
    import scipy

    file_path = images.path
    if isinstance(file_path, list):
        file_path = file_path[0]

    images = images.data

    # np.arrayをmmapへ変換
    order = 'C'
    dims = images.shape[1:]
    T = images.shape[0]
    shape_mov = (np.prod(dims), T)

    dir_path = join_file_path(file_path.split("/")[:-1])
    basename = file_path.split("/")[-1]
    fname_tot = memmap_frames_filename(basename, dims, T, order)

    mmap_images = np.memmap(
        join_file_path([dir_path, fname_tot]),
        mode='w+',
        dtype=np.float32,
        shape=prepare_shape(shape_mov),
        order=order)

    mmap_images = np.reshape(mmap_images.T, [T] + list(dims), order='F')
    mmap_images[:] = images[:]

    del images
    gc.collect()

    if params is None:
        ops = CNMFParams()
    else:
        ops = CNMFParams(params_dict=params)

    if 'dview' in locals():
        cm.stop_server(dview=dview)

    c, dview, n_processes = setup_cluster(
        backend='local', n_processes=None, single_thread=False)

    cnm = cnmf.CNMF(n_processes=n_processes, dview=dview, Ain=None, params=ops)
    cnm = cnm.fit(mmap_images)

    stop_server(dview=dview)

    # contours plot
    Cn = local_correlations(mmap_images.transpose(1, 2, 0))
    Cn[np.isnan(Cn)] = 0
    cnm.estimates.plot_contours(img=Cn)

    del mmap_images
    gc.collect()

    thr = params['thr']
    thr_method = 'nrg'
    swap_dim = False

    iscell = np.concatenate([
        np.ones(cnm.estimates.A.shape[-1]),
        np.zeros(cnm.estimates.b.shape[-1])
    ])

    cell_roi = get_roi(cnm.estimates.A, thr, thr_method, swap_dim, dims)
    cell_roi = np.stack(cell_roi)
    cell_roi = np.nanmax(cell_roi, axis=0).astype(float)
    cell_roi[cell_roi == 0] = np.nan

    non_cell_roi = get_roi(scipy.sparse.csc_matrix(cnm.estimates.b), thr, thr_method, swap_dim, dims)
    non_cell_roi = np.stack(non_cell_roi)
    non_cell_roi = np.nanmax(non_cell_roi, axis=0).astype(float)
    non_cell_roi[non_cell_roi == 0] = np.nan

    all_roi = np.nanmax(np.stack([cell_roi, non_cell_roi]), axis=0)

    # NWBの追加
    if nwbfile is not None:
        ### NWBにROIを追加
        roi_list = []
        n_cells = cnm.estimates.A.shape[-1]
        for i in range(n_cells):
            kargs = {}
            kargs['image_mask'] = cnm.estimates.A.T[i].T.toarray().reshape(dims)
            if hasattr(cnm.estimates, 'accepted_list'):
                kargs['accepted'] = i in cnm.estimates.accepted_list
            if hasattr(cnm.estimates, 'rejected_list'):
                kargs['rejected'] = i in cnm.estimates.rejected_list
            roi_list.append(kargs)

        ### backgroundsを追加
        bg_list = []
        for bg in cnm.estimates.b.T:
            kargs = {}
            kargs['image_mask'] = bg.reshape(dims)
            if hasattr(cnm.estimates, 'accepted_list'):
                kargs['accepted'] = False
            if hasattr(cnm.estimates, 'rejected_list'):
                kargs['rejected'] = False
            bg_list.append(kargs)

        nwbfile[NWBDATASET.ROI] = {
            'roi_list': roi_list,
            'bg_list': bg_list,
        }

        ### iscellを追加
        nwbfile[NWBDATASET.COLUMN] = {
            'roi_column': {
                'name': 'iscell',
                'discription': 'two columns - iscell & probcell',
                'data': iscell,
            }
        }

        ### Fluorescence
        n_rois = cnm.estimates.A.shape[-1]
        n_bg = len(cnm.estimates.f)

        nwbfile[NWBDATASET.FLUORESCENCE] = {
            'RoiResponseSeries': {
                'table_name': 'ROIs',
                'region': list(range(n_rois)),
                'name': 'RoiResponseSeries',
                'data': cnm.estimates.C.T,
                'unit': 'lumens',
            },
            'Background_Fluorescence_Response': {
                'table_name': 'Background',
                'region': list(range(n_rois, n_rois+n_bg)),
                'name': 'Background_Fluorescence_Response',
                'data': cnm.estimates.f.T,
                'unit': 'lumens',
            }
        }

    fluorescence = np.concatenate([
        cnm.estimates.C,
        cnm.estimates.f,
    ])

    info = {}
    info['images'] = ImageData(np.array(Cn * 255, dtype=np.uint8), file_name='images')
    info['fluorescence'] = FluoData(fluorescence, file_name='fluorescence')
    info['iscell'] = IscellData(iscell, file_name='iscell')
    info['cell_roi'] = RoiData(cell_roi, file_name='cell_roi')
    info['non_cell_roi'] = RoiData(non_cell_roi, file_name='non_cell_roi')
    info['all_roi'] = RoiData(all_roi, file_name='all_roi')
    info['nwbfile'] = nwbfile

    return info


if __name__ == '__main__':
    import os
    import numpy as np
    from motion_correction import caiman_mc
    import caiman

    info = {}
    file_path = join_file_path([
        '/Users', 'shogoakiyama', 'caiman_data', 
        'example_movies', 'Sue_2x_3000_40_-46.tif'])
    info['caiman_mc'] = caiman_mc(file_path)
    info['caiman_cnmf'] = caiman_cnmf(info['caiman_mc']['images'])
