from optinist.api.dataclass.dataclass import *
from optinist.wrappers.caiman_wrapper.cnmf import get_roi



def excute_delete_roi(node_dirpath, ids):
    import numpy as np
    
    from .utils import get_roi_manual
    from scipy.sparse import csc_matrix

    estimates = np.load(f'{node_dirpath}/caiman_cnmf.npy', allow_pickle=True).item()

    dims = estimates['dims']
    A = estimates.get('A').toarray()
    A_manual = estimates.get('A_manual', np.empty((dims[0]*dims[1], 0)))

    # get deleting ROI IDs
    deleting_ROI_IDs = []
    deleting_ROIs_manual_ID = []
    [deleting_ROI_IDs.append(id) for id in ids if id < A.shape[-1]]
    [deleting_ROIs_manual_ID.append(id - A.shape[-1]) for id in ids if id >= A.shape[-1]]
    
    # delete deleting ROIs from A and A_manual
    A = np.delete(A, deleting_ROI_IDs, axis=1)
    A_manual = np.delete(A_manual, deleting_ROIs_manual_ID, axis=1)

    # delete fluorescence
    estimates['C'] = np.delete(estimates['C'], ids, axis=0)

    thr = 0.9
    thr_method = 'nrg'
    swap_dim = False

    A = csc_matrix(A)
    cell_roi = get_roi(A, thr, thr_method, swap_dim, dims)
    start_idx = A.shape[1]
    if A_manual.shape[1] > 0:
        cell_roi = cell_roi + (get_roi_manual(A_manual, dims, start_idx))
    if len(cell_roi) == 0:
        cell_roi = np.zeros(dims)
    else:
        cell_roi = np.stack(cell_roi)
        cell_roi = np.nanmax(cell_roi, axis=0).astype(float)
    cell_roi[cell_roi == 0] = np.nan

    fluorescence = np.concatenate([
        estimates['C'],
        estimates['f'],
    ])

    estimates['A'] = A
    estimates['A_manual'] = A_manual

    info = {
        'fluorescence': FluoData(fluorescence, file_name='fluorescence'),
        'cell_roi': RoiData(cell_roi, file_name='cell_roi'),
        'cnmf_data': CaimanCnmfData(estimates), # save estimates object for further analysis
    }

    for v in info.values():
        if isinstance(v, BaseData):
            v.save_json(node_dirpath)