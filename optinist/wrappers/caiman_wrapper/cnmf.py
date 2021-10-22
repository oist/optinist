def caiman_cnmf(info, opts=None):
    images = info['images']
    from caiman import local_correlations
    from caiman.source_extraction.cnmf import cnmf
    from caiman.source_extraction.cnmf.params import CNMFParams
    import caiman.utils.visualization as visualization
    import numpy as np

    info = {}

    if opts is None:
        opts = CNMFParams()
    else:
        opts.change_params(params_dict=opts)

    cnm = cnmf.CNMF(1, params=opts, dview=None)
    cnm = cnm.fit(images)

    # contours plot
    Cn = local_correlations(images.transpose(1, 2, 0))
    Cn[np.isnan(Cn)] = 0
    cnm.estimates.plot_contours(img=Cn)

    # get roi center
    cont = visualization.get_contours(cnm.estimates.A, cnm.dims, thr=0.9, thr_method='nrg', swap_dim=False)
    cont_cent = np.zeros([len(cont), 2])
    for i in range(len(cont)):
        cont_cent[i, :] = np.nanmean(cont[i]['coordinates'], axis=0)

    iscell = np.zeros(cont_cent.shape[0])
    iscell[cnm.estimates.idx_components]=1

    info['images'] = np.array(Cn * 255, dtype=np.uint8)
    info['fluo'] = cnm.estimates.C
    info['iscell'] = iscell
    info['roi'] = cont_cent

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
