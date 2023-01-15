import numpy as np
import os
import datetime
import time
# from suite2p.gui.drawroi import masks_and_traces
from optinist.api.dataclass.dataclass import *
from scipy import stats
from suite2p import extraction, detection, classification, ROI

def masks_and_traces(ops, stat_manual, stat_orig):
    ''' main extraction function
        inputs: ops and stat
        creates cell and neuropil masks and extracts traces
        returns: F (ROIs x time), Fneu (ROIs x time), F_chan2, Fneu_chan2, ops, stat
        F_chan2 and Fneu_chan2 will be empty if no second channel
    '''
    if 'aspect' in ops:
        dy, dx = int(ops['aspect'] * 10), 10
    else:
        d0 = ops['diameter']
        dy, dx = (d0, d0) if isinstance(d0, int) else d0
    t0 = time.time()
    # Concatenate stat so a good neuropil function can be formed
    stat_all = stat_manual.copy()
    for n in range(len(stat_orig)):
        stat_all.append(stat_orig[n])
    stat_all = detection.stats.roi_stats(stat_all, dy, dx, ops['Ly'], ops['Lx'])
    cell_masks = [
        extraction.masks.create_cell_mask(stat, Ly=ops['Ly'], Lx=ops['Lx'], allow_overlap=ops['allow_overlap']) for stat in stat_all
    ]
    cell_pix = extraction.masks.create_cell_pix(stat_all, Ly=ops['Ly'], Lx=ops['Lx'])
    manual_roi_stats = stat_all[:len(stat_manual)]
    manual_cell_masks = cell_masks[:len(stat_manual)]
    manual_neuropil_masks = extraction.masks.create_neuropil_masks(
        ypixs=[stat['ypix'] for stat in manual_roi_stats],
        xpixs=[stat['xpix'] for stat in manual_roi_stats],
        cell_pix=cell_pix,
        inner_neuropil_radius=ops['inner_neuropil_radius'],
        min_neuropil_pixels=ops['min_neuropil_pixels'],
    )
    print('Masks made in %0.2f sec.' % (time.time() - t0))

    F, Fneu, F_chan2, Fneu_chan2, ops =extraction.extract_traces_from_masks(ops, 
                                                                  manual_cell_masks, 
                                                                  manual_neuropil_masks)

    # compute activity statistics for classifier
    npix = np.array([stat_orig[n]['npix'] for n in range(len(stat_orig))]).astype('float32')
    for n in range(len(manual_roi_stats)):
        manual_roi_stats[n]['npix_norm'] = manual_roi_stats[n]['npix'] / np.mean(npix[:100])  # What if there are less than 100 cells?
        manual_roi_stats[n]['compact'] = 1
        manual_roi_stats[n]['footprint'] = 2
        manual_roi_stats[n]['manual'] = 1  # Add manual key

    # subtract neuropil and compute skew, std from F
    dF = F - ops['neucoeff'] * Fneu
    sk = stats.skew(dF, axis=1)
    sd = np.std(dF, axis=1)

    for n in range(F.shape[0]):
        manual_roi_stats[n]['skew'] = sk[n]
        manual_roi_stats[n]['std'] = sd[n]
        manual_roi_stats[n]['med'] = [np.mean(manual_roi_stats[n]['ypix']), np.mean(manual_roi_stats[n]['xpix'])]

    dF = extraction.preprocess(
                F=dF,
                baseline=ops['baseline'],
                win_baseline=ops['win_baseline'],
                sig_baseline=ops['sig_baseline'],
                fs=ops['fs'],
                prctile_baseline=ops['prctile_baseline']
            )
    spks = extraction.dcnv.oasis(F=dF, batch_size=ops['batch_size'], tau=ops['tau'], fs=ops['fs'])

    return F, Fneu, F_chan2, Fneu_chan2, spks, ops, manual_roi_stats

def add_ROI(ops:dict, pos: list):
    def position(ops, pos):
        
        posx, posy, sizex, sizey = pos
        xrange = (np.arange(-1 * int(sizex), 1) + int(posx)).astype(np.int32)
        yrange = (np.arange(-1 * int(sizey), 1) + int(posy)).astype(np.int32)
        yrange += int(np.floor(sizey / 2)) + 1

        ellipse = np.zeros((yrange.size, xrange.size), np.bool)
        x, y = np.meshgrid(np.arange(0, xrange.size, 1), np.arange(0, yrange.size, 1))
        ellipse = ((y - sizey / 2) ** 2 / (sizey / 2) ** 2+ (x - sizex / 2) ** 2 / (sizex / 2) ** 2) <= 1

        ellipse = ellipse[:, np.logical_and(xrange >= 0, xrange < ops['Lx'])]
        xrange = xrange[np.logical_and(xrange >= 0, xrange < ops['Lx'])]
        ellipse = ellipse[np.logical_and(yrange >= 0, yrange < ops['Ly']), :]
        yrange = yrange[np.logical_and(yrange >= 0, yrange < ops['Ly'])]

        return ellipse, xrange, yrange

    ellipse, xrange, yrange = position(ops, pos)

    ops, stat, F, _, _, _ = extraction.create_masks_and_extract(ops, stat=ops['stat'])


    med = pos[2:4]
    x, y = np.meshgrid(xrange, yrange)
    ypix = y[ellipse].flatten()
    xpix = x[ellipse].flatten()
    lam = np.ones(ypix.shape)

    stat0 = [{'ypix': ypix, 'xpix': xpix, 'lam': lam, 'npix': ypix.size, 'med': med}]

    if not os.path.isfile(ops['reg_file']):
        ops['reg_file'] = os.path.join(ops.basename, 'data.bin')
    Fcell, Fneu, F_chan2, Fneu_chan2, spks, ops, stat_cell = masks_and_traces(ops, stat0, stat)


    stat_all = stat_cell.copy()
    F_all = F.copy()
    
    # iscell_prob = np.concatenate((iscell[:, np.newaxis], self.parent.probcell[:, np.newaxis]), axis=1)

    for n in range(len(stat)):
        stat_all.append(stat[n])
    F_all = np.concatenate((Fcell, F_all), axis=0)

    ops_classfile = ops.get('classifier_path')
    builtin_classfile = classification.builtin_classfile
    user_classfile = classification.user_classfile
    if ops_classfile:
        print(f'NOTE: applying classifier {str(ops_classfile)}')
        classfile = ops_classfile
    elif ops['use_builtin_classifier'] or not user_classfile.is_file():
        print(f'NOTE: Applying builtin classifier at {str(builtin_classfile)}')
        classfile = builtin_classfile
    else:
        print(f'NOTE: applying default {str(user_classfile)}')
        classfile = user_classfile

    iscell = classification.classify(stat=stat_all, classfile=classfile)
    iscell = iscell[:, 0].astype(bool)
    iscell[0] = True

    arrays = []
    for i, s in enumerate(stat_all):
        array = ROI(ypix=s['ypix'], xpix=s['xpix'], lam=s['lam'], med=s['med'], do_crop=False
        ).to_array(Ly=ops['Ly'], Lx=ops['Lx'])
        array *= i + 1
        arrays.append(array)
    im = np.stack(arrays)
    im[im == 0] = np.nan

    ops['stat'] = stat_all
    ops['F'] = F_all
    # ops['Fneu'] = Fneu

    info = {
        'ops': Suite2pData(ops),
        # 'max_proj': ImageData(ops['max_proj'], file_name='max_proj'),
        # 'Vcorr': ImageData(ops['Vcorr'], file_name='Vcorr'),
        'fluorescence': FluoData(F_all, file_name='fluorescence'),
        # 'iscell': IscellData(iscell, file_name='iscell'),
        'all_roi': RoiData(np.nanmax(im, axis=0), file_name='all_roi'), 
        # 'non_cell_roi': RoiData(np.nanmax(im[~iscell], axis=0), file_name='noncell_roi'),
        'cell_roi': RoiData(np.nanmax(im[iscell], axis=0), file_name='cell_roi'),
        # 'nwbfile': nwbfile,
    }

    return info
