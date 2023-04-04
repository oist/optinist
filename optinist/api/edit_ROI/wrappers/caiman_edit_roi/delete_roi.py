from optinist.api.dataclass.dataclass import *


def excute_delete_roi(node_dirpath, ids):
    import numpy as np
    from scipy.sparse import csc_matrix

    from .utils import get_roi

    # load data
    estimates = np.load(f'{node_dirpath}/caiman_cnmf.npy', allow_pickle=True).item()
    dims = estimates['dims']
    A = estimates.get('A').toarray()
    ims = estimates.get('ims')
    images = estimates.get('images')

    # delete ROI
    A = np.delete(A, ids, axis=1)
    ims = np.delete(ims, ids, axis=0)
    
    # delete fluorescence
    estimates['C'] = np.delete(estimates['C'], ids, axis=0)

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