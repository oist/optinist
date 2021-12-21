from wrappers.data_wrapper import *
from wrappers.args_check import args_check


@args_check
def caiman_cnmf(images: ImageData, params: dict=None) -> {'fluo': TimeSeriesData, 'iscell': IscellData, 'roi': RoiData}:
    import caiman
    from caiman import local_correlations
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
    fname_tot = caiman.paths.memmap_frames_filename(basename, dims, T, order)
    mmap_images = np.memmap(
        os.path.join(dir_path, fname_tot),
        mode='w+',
        dtype=np.float32,
        shape=caiman.mmapping.prepare_shape(shape_mov),
        order=order)

    mmap_images = np.reshape(mmap_images.T, [T] + list(dims), order='F')
    mmap_images[:] = images[:]

    if params is None:
        params = CNMFParams()
    else:
        params = CNMFParams(params_dict=params)

    cnm = cnmf.CNMF(1, params=params, dview=None)
    cnm = cnm.fit(mmap_images)

    # contours plot
    Cn = local_correlations(mmap_images.transpose(1, 2, 0))
    Cn[np.isnan(Cn)] = 0
    cnm.estimates.plot_contours(img=Cn)

    # get roi center
    cont = visualization.get_contours(cnm.estimates.A, cnm.dims, thr=0.9, thr_method='nrg', swap_dim=False)
    cont_cent = np.zeros([len(cont), 2])
    for i in range(len(cont)):
        cont_cent[i, :] = np.nanmean(cont[i]['coordinates'], axis=0)

    iscell = np.zeros(cont_cent.shape[0])
    iscell[cnm.estimates.idx_components] = 1

    info = {}
    info['images'] = ImageData(np.array(Cn * 255, dtype=np.uint8), func_name='caiman_cnmf', file_name='images')
    info['fluo'] = TimeSeriesData(cnm.estimates.C, func_name='caiman_cnmf', file_name='fluo')
    info['iscell'] = IscellData(iscell, func_name='caiman_cnmf', file_name='iscell')
    info['roi'] = RoiData(cont_cent)

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
