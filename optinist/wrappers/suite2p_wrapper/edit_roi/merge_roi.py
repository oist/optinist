import numpy as np
from optinist.api.dataclass.dataclass import *
from suite2p.detection.stats import median_pix, roi_stats
from .add_roi import get_im

def execute_merge_roi(node_dirpath, merged_roi_ids):
    ops = np.load(os.path.join(node_dirpath, 'suite2p.npy'), allow_pickle=True).item()
    iscell = np.load(os.path.join(node_dirpath, 'iscell.npy'))
    f_cell = ops['F']
    stat_orig=ops['stat']

    print('merging activity... this may take some time')
    ypix = np.zeros((0,), np.int32)
    xpix = np.zeros((0,), np.int32)
    lam = np.zeros((0,), np.float32)
    footprints  = np.array([])
    F    = np.zeros((0,f_cell.shape[1]), np.float32)
    # Fneu = np.zeros((0,f_cell.shape[1]), np.float32)
    merged_cells = []
    for n in np.array(merged_roi_ids):
        merged_cells.append(n)
    merged_cells = np.unique(np.array(merged_cells))

    for n in merged_cells:
        iscell[n] = False
        ypix = np.append(ypix, stat_orig[n]["ypix"])
        xpix = np.append(xpix, stat_orig[n]["xpix"])
        lam = np.append(lam, stat_orig[n]["lam"])
        footprints = np.append(footprints, stat_orig[n]["footprint"])
        F    = np.append(F, f_cell[n,:][np.newaxis,:], axis=0)
        # Fneu = np.append(Fneu, f_neu[n,:][np.newaxis,:], axis=0)

    # remove overlaps
    ipix = np.concatenate((ypix[:,np.newaxis], xpix[:,np.newaxis]), axis=1)
    _, goodi = np.unique(ipix, return_index=True, axis=0)
    ypix = ypix[goodi]
    xpix = xpix[goodi]
    lam = lam[goodi]

    ### compute statistics of merges
    stat0 = {}
    stat0['imerge'] = merged_cells
    if 'iplane' in stat_orig[merged_cells[0]]:
        stat0['iplane'] = stat_orig[merged_cells[0]]['iplane']
    stat0['ypix'] = ypix
    stat0['xpix'] = xpix
    stat0['med'] = median_pix(ypix, xpix)
    stat0['lam'] = lam / lam.sum()
    
    stat0['chan2_prob'] = -1
    stat0["inmerge"] = -1

    ### compute activity of merged cells
    F = F.mean(axis=0)
    # Fneu = Fneu.mean(axis=0)


    if 'aspect' in ops:
        d0 = np.array([int(ops['aspect']*10), 10])
    else:
        d0 = ops['diameter']
        if isinstance(d0, int):
            d0 = [d0, d0]

    # add cell to structs
    ly = ops['Ly']
    lx = ops['Lx']
    stat_orig = np.concatenate((stat_orig, np.array([stat0])), axis=0)
    iscell = np.append(iscell, [True])
    stat_orig = roi_stats(stat_orig, d0[0], d0[1], ly, lx)
    stat_orig[-1]['lam'] = stat_orig[-1]['lam'] * merged_cells.size
    f_cell = np.concatenate((f_cell, F[np.newaxis,:]), axis=0)
    # f_neu = np.concatenate((f_neu, Fneu[np.newaxis,:]), axis=0)
    ops['stat'] = stat_orig
    ops['F'] = f_cell
    im = get_im(ops, stat_orig)

    info = {
        'ops': Suite2pData(ops),
        'fluorescence': FluoData(f_cell, file_name='fluorescence'),
        'iscell': IscellData(iscell, file_name='iscell'),
        'all_roi': RoiData(np.nanmax(im, axis=0), file_name='all_roi'), 
        'non_cell_roi': RoiData(np.nanmax(im[~iscell], axis=0), file_name='noncell_roi'),
        'cell_roi': RoiData(np.nanmax(im[iscell], axis=0), file_name='cell_roi'),
    }
    
    for k, v in info.items():
        if isinstance(v, BaseData):
            v.save_json(node_dirpath)

    cell_roi_data = np.where(np.isnan(info['cell_roi'].data), None, info['cell_roi'].data)
    max_index = len(info['fluorescence'].data)
    return cell_roi_data.tolist(), max_index
