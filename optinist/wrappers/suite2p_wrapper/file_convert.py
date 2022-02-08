from wrappers.data_wrapper import *
from wrappers.args_check import args_check
import gc

from cui_api.const import BASE_DIR
from cui_api.utils import join_file_path


def suite2p_file_convert(
        image: ImageData, nwbfile: NWBFile=None, params: dict=None
    ) -> {'ops': Suite2pData, 'images': ImageData}:
    import os
    import numpy as np
    from natsort import natsorted
    from suite2p import io, default_ops
    print('start suite2_file_convert')

    file_path = image.path
    data_path = '/'.join(file_path.split('/')[:-1])
    data_name = file_path.split('/')[-1]
    print(file_path)
    print(data_path)
    print(data_name)
    ### data pathと保存pathを指定
    db = {
        'data_path': [data_path],
        'tiff_list': [data_name],
        'save_path0': BASE_DIR,
        'save_folder': 'suite2p'
    }

    ops = {**default_ops(), **params, **db}

    ops['input_format'] = 'tif'

    # save folderを指定
    save_folder = join_file_path([ops['save_path0'], ops['save_folder']])
    os.makedirs(save_folder, exist_ok=True)
    plane_folders = natsorted([
        f.path for f in os.scandir(save_folder) 
        if f.is_dir() and f.name[:5]=='plane'])

    ops_path = [join_file_path([f, 'ops.npy']) for f in plane_folders]

    # copy file format to a binary file
    convert_funs = {
        'h5': io.h5py_to_binary,
        'sbx': io.sbx_to_binary,
        'mesoscan': io.mesoscan_to_binary,
        'haus': lambda ops: haussio.load_haussio(ops['data_path'][0]).tosuite2p(ops.copy()),
        'bruker': io.ome_to_binary,
        'tif': io.tiff_to_binary
    }

    # save ops.npy(parameter) and data.bin
    ops = convert_funs[ops['input_format']](ops.copy())

    info = {}
    images = io.BinaryFile(read_filename=ops['reg_file'], Ly=ops['Ly'], Lx=ops['Lx']).data
    info['images'] = ImageData(images, func_name='suite2p_registration', file_name='images')
    info['meanImg'] = ImageData(ops['meanImg'], func_name='suite2p_convert', file_name='meanImg')
    info['ops'] = Suite2pData(ops, func_name='suite2p_convert', file_name='ops')

    del images, io
    gc.collect()

    return info
