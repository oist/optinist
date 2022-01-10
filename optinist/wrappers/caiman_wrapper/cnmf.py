from wrappers.data_wrapper import *
from wrappers.args_check import args_check


@args_check
def caiman_cnmf(
        images: ImageData, nwbfile: NWBFile=None, params: dict=None
    ) -> {'fluo': TimeSeriesData, 'iscell': IscellData, 'roi': RoiData}:
    import caiman
    from caiman import local_correlations, stop_server
    from caiman.paths import memmap_frames_filename
    from caiman.mmapping import prepare_shape
    from caiman.cluster import setup_cluster
    from caiman.source_extraction.cnmf import cnmf
    from caiman.source_extraction.cnmf.params import CNMFParams
    import caiman.utils.visualization as visualization
    import numpy as np

    file_path = images.path
    images = images.data

    # np.arrayをmmapへ変換
    order = 'C'
    dims = images.shape[1:]
    T = images.shape[0]
    shape_mov = (np.prod(dims), T)

    dir_path = os.path.dirname(file_path)
    basename = os.path.splitext(os.path.basename(file_path))[0]
    fname_tot = memmap_frames_filename(basename, dims, T, order)
    mmap_images = np.memmap(
        os.path.join(dir_path, fname_tot),
        mode='w+',
        dtype=np.float32,
        shape=prepare_shape(shape_mov),
        order=order)

    mmap_images = np.reshape(mmap_images.T, [T] + list(dims), order='F')
    mmap_images[:] = images[:]

    if params is None:
        params = CNMFParams()
    else:
        params = CNMFParams(params_dict=params)

    c, dview, n_processes = setup_cluster(
        backend='local', n_processes=None, single_thread=False)

    cnm = cnmf.CNMF(n_processes, params=params, dview=dview)
    cnm = cnm.fit(mmap_images)

    stop_server(dview=dview)

    # contours plot
    Cn = local_correlations(mmap_images.transpose(1, 2, 0))
    Cn[np.isnan(Cn)] = 0
    cnm.estimates.plot_contours(img=Cn)

    # get roi center
    cont = visualization.get_contours(
        cnm.estimates.A, cnm.dims, thr=0.9, thr_method='nrg', swap_dim=False)
    cont_cent = np.zeros([len(cont), 2])
    for i in range(len(cont)):
        cont_cent[i, :] = np.nanmean(cont[i]['coordinates'], axis=0)

    iscell = np.zeros(cont_cent.shape[0])
    iscell[cnm.estimates.idx_components] = 1

    # NWBの追加

    ### NWBにROIを追加
    roi_list = []
    n_cells = cnm.estimates.A.shape[-1]
    for i in range(n_cells):
        kargs = {}
        kargs['image_mask'] = cnm.estimates.A.T[i].T.toarray().reshape(cnm.estimates.dims)
        # image_mask_list.append()
        if hasattr(cnm.estimates, 'accepted_list'):
            kargs['accepted'] = i in cnm.estimates.accepted_list
        if hasattr(cnm.estimates, 'rejected_list'):
            kargs['rejected'] = i in cnm.estimates.rejected_list
        roi_list.append(kargs)

    nwbfile = nwb_add_roi(nwbfile, roi_list)

    ### backgroundsを追加
    bg_list = []
    for bg in cnm.estimates.b.T:
        kwargs = dict(
            image_mask=bg.reshape(cnm.estimates.dims),
            accepted=False,
            rejected=False
        )
        bg_list.append(kargs)

    nwbfile = nwb_add_roi(nwbfile, bg_list)

    ### iscellを追加
    nwbfile = nwb_add_column(
        nwbfile, 'iscell', 'two columns - iscell & probcell', iscell)

    ### Fluorescence
    starting_time = 0
    imaging_rate = nwbfile.imaging_planes['ImagingPlane'].imaging_rate
    timestamps = np.arange(cnm.estimates.f.shape[1]) / imaging_rate + starting_time
    n_rois = cnm.estimates.A.shape[-1]
    n_bg = len(cnm.estimates.f)

    ### roiを追加
    nwbfile = nwb_add_fluorescence(
        nwbfile,
        table_name='ROIs',
        region=list(range(n_rois)),
        name='RoiResponseSeries',
        data=cnm.estimates.C.T,
        unit='lumens',
        timestamps=timestamps
    )

    ### backgroundsを追加
    nwbfile = nwb_add_fluorescence(
        nwbfile,
        table_name='Background',
        region=list(range(n_rois, n_rois+n_bg)),
        name='Background_Fluorescence_Response',
        data=cnm.estimates.f.T,
        unit='lumens',
        timestamps=timestamps
    )

    info = {}
    info['images'] = ImageData(np.array(Cn * 255, dtype=np.uint8), func_name='caiman_cnmf', file_name='images')
    info['F'] = TimeSeriesData(cnm.estimates.C, func_name='caiman_cnmf', file_name='F')
    info['iscell'] = IscellData(iscell, func_name='caiman_cnmf', file_name='iscell')
    info['roi'] = RoiData(cont_cent)
    info['nwbfile'] = nwbfile

    return info


if __name__ == '__main__':
    import os
    import numpy as np
    from motion_correction import caiman_mc
    import caiman

    info = {}
    file_path = os.path.join(
        '/Users', 'shogoakiyama', 'caiman_data', 
        'example_movies', 'Sue_2x_3000_40_-46.tif')
    info['caiman_mc'] = caiman_mc(file_path)
    info['caiman_cnmf'] = caiman_cnmf(info['caiman_mc']['images'])
