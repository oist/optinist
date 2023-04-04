from optinist.api.dataclass.dataclass import *


def execute_add_ROI(node_dirpath, posx, posy, sizex, sizey):
    import numpy as np
    from scipy.sparse import csc_matrix

    from .utils import create_mask, get_roi

    # load data
    estimates = np.load(f'{node_dirpath}/caiman_cnmf.npy', allow_pickle=True).item()
    dims = estimates['dims']
    A = estimates.get('A').toarray()
    ims = estimates.get('ims')
    images = estimates.get('images')

    # Create mask for new roi and append to A, ims
    new_roi = create_mask(posx, posy, sizex, sizey, dims)
    A = np.hstack((new_roi.reshape(-1, 1), A))
    ims = np.concatenate((new_roi[np.newaxis, :, :], ims), axis=0)

    # extract fluorescence of new ROI
    num_frames = images.shape[0]
    reshapeImages = images.reshape([images.shape[2] * images.shape[1], num_frames])
    new_fluorescence = np.zeros(num_frames)
    new_fluorescence = np.mean(reshapeImages[new_roi.reshape(-1, 1)[:,0] > 0], axis=0)
    estimates['C'] = np.vstack([new_fluorescence, estimates['C']])

    cell_roi = get_roi(ims)
    fluorescence = np.concatenate([
        estimates['C'],
        estimates['f'],
    ])

    # save data
    estimates['A'] = csc_matrix(A)
    estimates['ims'] = ims
    # breakpoint()
    info = {
        'fluorescence': FluoData(fluorescence, file_name='fluorescence'),
        'cell_roi': RoiData(cell_roi, file_name='cell_roi'),
        'cnmf_data': CaimanCnmfData(estimates), # save estimates object for further analysis
    }

    for v in info.values():
        if isinstance(v, BaseData):
            v.save_json(node_dirpath)