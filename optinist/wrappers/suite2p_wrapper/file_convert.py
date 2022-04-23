from optinist.api.dataclass.dataclass import *

from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import join_filepath


def suite2p_file_convert(
        image: ImageData,
        params: dict=None
    ) -> dict(ops=Suite2pData):
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
        'save_path0': DIRPATH.OPTINIST_DIR,
        'save_folder': 'suite2p'
    }

    ops = {**default_ops(), **params, **db}

    # save folderを指定
    create_directory(
        join_filepath([ops['save_path0'], ops['save_folder']])
    )

    # save ops.npy(parameter) and data.bin
    ops = io.tiff_to_binary(ops.copy())

    info = {
        'meanImg': ImageData(ops['meanImg'], file_name='meanImg'),
        'ops': Suite2pData(ops, file_name='ops')
    }

    return info
