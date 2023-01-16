import numpy as np
from optinist.api.dataclass.dataclass import *
from .utils import get_im, save_json_data

def excute_delete_roi(node_dirpath, delete_roi_ids):
    ops = np.load(os.path.join(node_dirpath, 'suite2p.npy'), allow_pickle=True).item()
    iscell = ops.get('iscell')
    stat = ops.get('stat')

    for id in delete_roi_ids:
        iscell[id] = False

    im = get_im(ops, stat)

    ops['iscell'] = iscell

    info = save_json_data(ops, im, save_path=node_dirpath, save_data=['ops', 'non_cell_roi', 'cell_roi'])    

    cell_roi_data = info['cell_roi'].data
    cell_roi_data = np.where(np.isnan(cell_roi_data), None, cell_roi_data)
    max_index = len(ops['F'])
    return cell_roi_data.tolist(), max_index