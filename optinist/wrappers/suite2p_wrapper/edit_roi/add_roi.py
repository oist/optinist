import numpy as np
import os
from optinist.api.dataclass.dataclass import *
from .utils import get_im, get_stat0_add_roi, masks_and_traces, save_json_data

def execute_add_ROI(node_dirpath, pos: list):
    ops = np.load(os.path.join(node_dirpath, 'suite2p.npy'), allow_pickle=True).item()
    iscell = ops.get('iscell')
    stat = ops.get('stat')
    stat0 = get_stat0_add_roi(ops, pos)

    Fcell, Fneu, F_chan2, Fneu_chan2, spks, ops, stat_cell = masks_and_traces(ops, stat0, ops['stat'])
    
    stat_all = stat_cell.copy()
    for n in range(len(stat)):
        stat_all.append(stat[n])

    # Save fluorescence traces
    F_all = np.concatenate((Fcell, ops['F']), axis=0)
    Fneu_all = np.concatenate((Fneu, ops['Fneu']), axis=0)

    iscell = np.append(True, iscell)
    im = get_im(ops, stat_all)

    ops['stat'] = stat_all
    ops['F'] = F_all
    ops['Fneu'] = Fneu_all
    ops['iscell'] = iscell

    info = save_json_data(ops, im, save_path=node_dirpath, save_data=['ops', 'fluorescence', 'all_roi', 'non_cell_roi', 'cell_roi'])

    max_index = len(info['fluorescence'].data)
    return max_index

