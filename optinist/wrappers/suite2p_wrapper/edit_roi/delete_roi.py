import numpy as np
from optinist.api.dataclass.dataclass import *
from .add_roi import get_im

def excute_delete_roi(node_dirpath, delete_roi_ids):
    ops = np.load(os.path.join(node_dirpath, 'suite2p.npy'), allow_pickle=True).item()
    iscell = np.load(os.path.join(node_dirpath, 'iscell.npy'))
    stat=ops['stat']

    for id in delete_roi_ids:
        iscell[id] = False

    im = get_im(ops, stat)
    info = {
        'iscell': IscellData(iscell, file_name='iscell'),
        'non_cell_roi': RoiData(np.nanmax(im[~iscell], axis=0), file_name='noncell_roi'),
        'cell_roi': RoiData(np.nanmax(im[iscell], axis=0), file_name='cell_roi'),
    }
    
    for k, v in info.items():
        if isinstance(v, BaseData):
            v.save_json(node_dirpath)

    cell_roi_data = info['cell_roi'].data
    cell_roi_data = np.where(np.isnan(cell_roi_data), None, cell_roi_data)
    max_index = len(ops['F'])
    return cell_roi_data.tolist(), max_index