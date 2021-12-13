from wrappers.data_wrapper import *
from wrappers.args_check import args_check


@args_check
def suite2p_roi(ops: Suite2pData, params: dict=None):
    ops = ops.data

    import numpy as np
    from suite2p import extraction, classification, detection, ROI
    from suite2p import default_ops

    ops = {**ops, **params, **default_ops()}

    # ROI detection
    classfile = classification.user_classfile
    ops, stat = detection.detect(ops=ops, classfile=classfile)

    ######## ROI EXTRACTION ##############
    ops, stat, F, Fneu, F_chan2, Fneu_chan2 = extraction.create_masks_and_extract(ops, stat)

    ######## ROI CLASSIFICATION ##############
    iscell = classification.classify(stat=stat, classfile=classfile)
    iscell = iscell[:, 0].astype(bool)

    arrays = []
    for i, s in enumerate(stat):
        array = ROI(
            ypix=s['ypix'], xpix=s['xpix'], lam=s['lam'], med=s['med'], do_crop=False
        ).to_array(Ly=ops['Ly'], Lx=ops['Lx'])
        array *= i + 1
        arrays.append(array)

    im = np.stack(arrays)
    im[im == 0] = np.nan

    ops['ROI_found'] = np.nanmax(im, axis=0)
    ops['non_cell_roi'] = np.nanmax(im[~iscell])
    ops['cell_roi'] = np.nanmax(im[iscell], axis=0)
    ops['F'] = F
    ops['Fneu'] = Fneu

    info = {}
    info['max_proj'] = ImageData(
        ops['max_proj'],
        func_name='suite2p_roi',
        file_name='max_proj'
    )
    import pdb; pdb.set_trace()
    info['F'] = TimeSeriesData(
        F,
        func_name='suite2p_roi',
        file_name='F'
    )

    info['ops'] = Suite2pData(ops)

    return info
