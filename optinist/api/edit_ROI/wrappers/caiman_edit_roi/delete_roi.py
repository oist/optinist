from optinist.api.dataclass.dataclass import *


def excute_delete_roi(node_dirpath, ids):
    import numpy as np

    from .utils import get_roi

    # load data
    cnmf_data = np.load(f'{node_dirpath}/caiman_cnmf.npy', allow_pickle=True).item()
    is_cell = cnmf_data.get('is_cell')
    ims = cnmf_data.get('ims')

    # delete ROI
    is_cell[ids] = False

    cell_roi = get_roi(ims)

    cnmf_data['ims'] = ims
    cnmf_data['is_cell'] = is_cell

    info = {
        'cell_roi': RoiData(np.nanmax(cell_roi[is_cell], axis=0), file_name='cell_roi',),
        'cnmf_data': CaimanCnmfData(cnmf_data),
    }

    for v in info.values():
        if isinstance(v, BaseData):
            v.save_json(node_dirpath)
