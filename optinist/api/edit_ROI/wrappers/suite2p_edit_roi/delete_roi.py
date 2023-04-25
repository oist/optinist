import numpy as np

from optinist.api.dataclass.dataclass import *
from optinist.api.edit_ROI.utils import save_edit_ROI_data

from .utils import set_nwbfile

@save_edit_ROI_data
def excute_delete_roi(node_dirpath, ids):
    ops = np.load(os.path.join(node_dirpath, 'suite2p.npy'), allow_pickle=True).item()
    iscell = ops.get('iscell')
    delete_roi = ops.get('delete_roi', [])
    im = ops.get('im')

    delete_roi = ops.get('delete_roi') if ops.get('delete_roi') else []
    for id in ids:
        iscell[id] = False
        delete_roi.append(id + 1)

    ops['iscell'] = iscell
    ops['delete_roi'] = delete_roi

    info = {
        'ops': Suite2pData(ops),
        'non_cell_roi': RoiData(
            np.nanmax(im[~iscell], axis=0), file_name='noncell_roi'
        ),
        'cell_roi': RoiData(np.nanmax(im[iscell], axis=0), file_name='cell_roi'),
        'nwbfile': set_nwbfile(ops),
    }
    
    return info
