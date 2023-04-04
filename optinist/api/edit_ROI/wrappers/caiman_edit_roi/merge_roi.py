from optinist.api.dataclass.dataclass import *


def execute_merge_roi(node_dirpath: str, ids: list):
    import numpy as np
    from scipy.sparse import csc_matrix

    from .utils import get_roi

    # load data
    estimates = np.load(f'{node_dirpath}/caiman_cnmf.npy', allow_pickle=True).item()
    dims = estimates['dims']
    A = estimates.get('A').toarray()
    ims = estimates.get('ims')
    images = estimates.get('images')

    # get merging ROI
    merging_ROIs = []
    [merging_ROIs.append(ims[id,:,:]) for id in ids]

    # delete merging ROIs from A and ims
    A = np.delete(A, ids, axis=1)
    ims = np.delete(ims, ids, axis=0)
    
    # get merged ROI and append to A and ims
    merged_ROI = np.maximum.reduce(merging_ROIs)
    A = np.hstack((A, merged_ROI.reshape(-1, 1)))
    ims = np.concatenate((ims, merged_ROI[np.newaxis, :, :]), axis=0)

    # get merged F
    merged_f = np.mean((estimates['C'][ids, :]), axis=0)
    estimates['C'] = np.delete(estimates['C'], ids, axis=0)
    estimates['C'] = np.vstack([estimates['C'], merged_f])

    cell_roi = get_roi(ims)

    fluorescence = np.concatenate([
        estimates['C'],
        estimates['f'],
    ])

    estimates['A'] = csc_matrix(A)
    estimates['ims'] = ims

    info = {
        'fluorescence': FluoData(fluorescence, file_name='fluorescence'),
        'cell_roi': RoiData(cell_roi, file_name='cell_roi'),
        'cnmf_data': CaimanCnmfData(estimates), # save estimates object for further analysis
    }

    for v in info.values():
        if isinstance(v, BaseData):
            v.save_json(node_dirpath)