from optinist.wrappers.data_wrapper import *

from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import join_filepath


def suite2p_file_convert(
        image: ImageData,
        nwbfile: NWBFile=None,
        params: dict=None
    ) -> dict(ops=Suite2pData):
    import os
    from suite2p import io, default_ops
    print('start suite2_file_convert')

    data_path_list = []
    data_name_list = []
    for file_path in image.path:
        data_path_list.append('/'.join(file_path.split('/')[:-1]))
        data_name_list.append(file_path.split('/')[-1])

    print(data_path_list)
    print(data_name_list)
    ### data pathと保存pathを指定
    db = {
        'data_path': data_path_list,
        'tiff_list': data_name_list,
        'save_path0': DIRPATH.BASE_DIR,
        'save_folder': 'suite2p'
    }

    ops = {**default_ops(), **params, **db}

    ops['input_format'] = 'tif'

    # save folderを指定
    save_folder = join_filepath([ops['save_path0'], ops['save_folder']])
    os.makedirs(save_folder, exist_ok=True)

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
    info['meanImg'] = ImageData(ops['meanImg'], file_name='meanImg')
    info['ops'] = Suite2pData(ops, file_name='ops')

    return info
