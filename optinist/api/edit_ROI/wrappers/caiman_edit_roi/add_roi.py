from optinist.api.dataclass.dataclass import *
from optinist.wrappers.caiman_wrapper.cnmf import get_roi

from .utils import get_roi_manual


def execute_add_ROI(node_dirpath, posx, posy, sizex, sizey):
    import numpy as np

    from .utils import create_mask

    estimates = np.load(f'{node_dirpath}/caiman_cnmf.npy', allow_pickle=True).item()

    dims = estimates['dims']
    A = estimates.get('A')
    A_manual = estimates.get('A_manual', np.empty((dims[0]*dims[1], 0)))

    new_roi = create_mask(posx, posy, sizex, sizey, dims).T
    A_manual = np.hstack((A_manual, new_roi.reshape(-1, 1)))

    images = estimates.get('images')
    num_frames = images.shape[0]
    reshapeImages = images.reshape([images.shape[2] * images.shape[1], num_frames])
    new_fluorescence = np.zeros(num_frames)
    new_fluorescence = np.mean(reshapeImages[new_roi.reshape(-1, 1)[:,0] > 0], axis=0)
    estimates['C'] = np.vstack([estimates['C'], new_fluorescence])

    thr = 0.9
    thr_method = 'nrg'
    swap_dim = False

    cell_roi = get_roi(A, thr, thr_method, swap_dim, dims)
    start_idx = A.shape[1]
    if A_manual.shape[1] > 0:
        cell_roi = cell_roi + (get_roi_manual(A_manual, dims, start_idx))
    cell_roi = np.stack(cell_roi)
    cell_roi = np.nanmax(cell_roi, axis=0).astype(float)
    cell_roi[cell_roi == 0] = np.nan

    fluorescence = np.concatenate([
        estimates['C'],
        estimates['f'],
    ])

    estimates['A_manual'] = A_manual

    info = {
        'fluorescence': FluoData(fluorescence, file_name='fluorescence'),
        'cell_roi': RoiData(cell_roi, file_name='cell_roi'),
        'cnmf_data': CaimanCnmfData(estimates), # save estimates object for further analysis
    }

    for v in info.values():
        if isinstance(v, BaseData):
            v.save_json(node_dirpath)