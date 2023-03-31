from optinist.api.dataclass.dataclass import *
from optinist.wrappers.caiman_wrapper.cnmf import get_roi


def execute_merge_roi(node_dirpath: str, ids: list):
    import numpy as np
    from scipy.sparse import csc_matrix

    from .utils import get_roi_manual

    estimates = np.load(f'{node_dirpath}/caiman_cnmf.npy', allow_pickle=True).item()

    dims = estimates['dims']
    A = estimates.get('A').toarray()
    A_manual = estimates.get('A_manual', np.empty((dims[0]*dims[1], 0)))

    # get merging ROI
    merging_ROIs = []
    merging_ROIs_manual = []
    [merging_ROIs.append(A[:, id]) for id in ids if id < A.shape[-1]]
    [merging_ROIs_manual.append(A_manual[:, id - A.shape[-1]]) for id in ids if id >= A.shape[-1]]

    # get merging ROI ID
    merging_ROI_IDs = []
    merging_ROIs_manual_ID = []
    [merging_ROI_IDs.append(id) for id in ids if id < A.shape[-1]]
    [merging_ROIs_manual_ID.append(id - A.shape[-1]) for id in ids if id >= A.shape[-1]]
    
    # delete merging ROIs from A and A_manual
    A = np.delete(A, merging_ROI_IDs, axis=1)
    A_manual = np.delete(A_manual, merging_ROIs_manual_ID, axis=1)

    # get merged ROI
    merged_ROI = np.maximum(*(merging_ROIs + merging_ROIs_manual))
    merged_ROI = merged_ROI.reshape(-1, 1)
    if not merging_ROIs_manual:
        A = np.hstack((A, merged_ROI))
    else:
        A_manual = np.hstack((A_manual, merged_ROI))

    # get merged F
    merged_f = np.mean((estimates['C'][ids, :]), axis=0)
    estimates['C'] = np.delete(estimates['C'], ids, axis=0)
    estimates['C'] = np.vstack([estimates['C'], merged_f])

    thr = 0.9
    thr_method = 'nrg'
    swap_dim = False

    A = csc_matrix(A)
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