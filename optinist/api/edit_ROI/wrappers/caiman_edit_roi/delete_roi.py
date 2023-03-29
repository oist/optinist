from optinist.api.dataclass.dataclass import *
from optinist.wrappers.caiman_wrapper.cnmf import get_roi



def excute_delete_roi(node_dirpath, ids):
    import numpy as np
    
    from scipy.sparse import csc_matrix

    estimates = np.load(f'{node_dirpath}/caiman_cnmf.npy', allow_pickle=True).item()

    dims = estimates['dims']
    A = estimates.get('A').toarray()

    A = np.delete(A, ids, axis=1)

    estimates['C'] = np.delete(estimates['C'], ids, axis=0)


    thr = 0.9
    thr_method = 'nrg'
    swap_dim = False

    iscell = np.concatenate([
        np.ones(A.shape[-1]),
        np.zeros(estimates['b'].shape[-1])
    ])
    A = csc_matrix(A)
    cell_roi = get_roi(A, thr, thr_method, swap_dim, dims)
    if len(cell_roi) == 0:
        cell_roi = np.zeros(dims)
    else:
        cell_roi = np.stack(cell_roi)
        cell_roi = np.nanmax(cell_roi, axis=0).astype(float)
    cell_roi[cell_roi == 0] = np.nan

    non_cell_roi = get_roi(csc_matrix(estimates['b']), thr, thr_method, swap_dim, dims)
    non_cell_roi = np.stack(non_cell_roi)
    non_cell_roi = np.nanmax(non_cell_roi, axis=0).astype(float)
    non_cell_roi[non_cell_roi == 0] = np.nan

    all_roi = np.nanmax(np.stack([cell_roi, non_cell_roi]), axis=0)

    fluorescence = np.concatenate([
        estimates['C'],
        estimates['f'],
    ])

    estimates['A'] = A

    info = {
        'fluorescence': FluoData(fluorescence, file_name='fluorescence'),
        'iscell': IscellData(iscell, file_name='iscell'),
        'all_roi': RoiData(all_roi, file_name='all_roi'),
        'cell_roi': RoiData(cell_roi, file_name='cell_roi'),
        'cnmf_data': CaimanCnmfData(estimates), # save estimates object for further analysis
    }

    for v in info.values():
        if isinstance(v, BaseData):
            v.save_json(node_dirpath)