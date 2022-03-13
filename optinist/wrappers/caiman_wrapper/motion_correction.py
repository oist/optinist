from wrappers.data_wrapper import *
from wrappers.args_check import args_check
import gc
from wrappers.nwb_wrapper.const import NWBDATASET


def caiman_mc(
        image: ImageData, nwbfile: NWBFile=None, params: dict=None
    ) -> {'mc_images': ImageData, 'iscell': IscellData}:
    file_path = image.path
    import numpy as np
    from caiman import load, save_memmap, load_memmap, stop_server
    from caiman.source_extraction.cnmf.params import CNMFParams
    from caiman.motion_correction import MotionCorrect
    from caiman.cluster import setup_cluster
    from caiman.base.rois import extract_binary_masks_from_structural_channel
    info = {}

    if params is None:
        opts = CNMFParams()
    else:
        opts = CNMFParams()
        opts.change_params(params_dict=params)

    c, dview, n_processes = setup_cluster(
        backend='local', n_processes=None, single_thread=False)

    mc = MotionCorrect(
        file_path, dview=dview, **opts.get_group('motion'))

    mc.motion_correct(save_movie=True)
    border_to_0 = 0 if mc.border_nan == 'copy' else mc.border_to_0

    # memory mapping
    fname_new = save_memmap(
        mc.mmap_file, base_name='memmap_', order='C',border_to_0=border_to_0)

    stop_server(dview=dview)

    # now load the file
    Yr, dims, T = load_memmap(fname_new)

    images = np.array(
        Yr.T.reshape((T,) + dims, order='F'))

    meanImg = images.mean(axis=0)
    iscell = extract_binary_masks_from_structural_channel(
        meanImg, gSig=7, expand_method='dilation')[0].reshape(
            meanImg.shape[0], meanImg.shape[1], -1).transpose(2, 0, 1)

    info['images'] = ImageData(images, func_name='caiman_mc', file_name='mc_images')
    info['meanImg'] = ImageData(meanImg, func_name='caiman_mc', file_name='meanImg')
    info['iscell'] = IscellData(meanImg, func_name='caiman_mc', file_name='iscell')

    xy_trans_data = (np.array(mc.x_shifts_els), np.array(mc.y_shifts_els)) \
                    if params['pw_rigid'] else np.array(mc.shifts_rig)

    if nwbfile is not None:
        nwbfile[NWBDATASET.MOTION_CORRECTION] = {
            'caiman_mc': {
                'mc_data': info['images'],
                'xy_trans_data': xy_trans_data,
            }
        }

    info['nwbfile'] = nwbfile

    return info


if __name__ == '__main__':
    import os
    from cui_api.utils import join_file_path
    info = {}
    file_path = join_file_path([
        '/Users', 'shogoakiyama', 'caiman_data', 
        'example_movies', 'Sue_2x_3000_40_-46.tif'])
    info['caiman_mc'] = caiman_mc(file_path)
    print(info)
