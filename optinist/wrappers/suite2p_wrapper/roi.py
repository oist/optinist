import numpy as np
from pathlib import Path
from suite2p import ROI


def roi(save_path, Ly, Lx):
    info = {}
    stats_file = Path(save_path).joinpath('stat.npy')
    iscell = np.load(
        Path(save_path).joinpath('iscell.npy'), allow_pickle=True
    )[:, 0].astype(bool)
    stats = np.load(stats_file, allow_pickle=True)

    label_id = True
    arrays = []
    for i, s in enumerate(stats):
        array = ROI(
            ypix=s['ypix'], xpix=s['xpix'], lam=s['lam'], med=s['med'], do_crop=False
        ).to_array(Ly=Ly, Lx=Lx)

        if label_id:
            array *= i + 1
        arrays.append(array)

    im = np.stack(arrays)
    im[im == 0] = np.nan

    info['im'] = im
    info['iscell'] = iscell

    return info


if __name__ == '__main__':
    import os
    from run_s2p import run_s2p
    file_path = os.path.join(
        '/Users', 'shogoakiyama', 'Desktop', 'optinist', 
        'optinist', 'data', 'Sue_2x_3000_40_-46.tif')
    info = run_s2p(file_path)
    roi(info['save_path'], info['Lx'], info['Ly'])
