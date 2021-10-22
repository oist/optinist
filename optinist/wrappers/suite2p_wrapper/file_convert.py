def suite2p_file_convert(ops, opts=None):
    import os
    import numpy as np
    from natsort import natsorted
    import suite2p.io as io
    from suite2p import default_ops

    path = ops['path']
    data_path = '/'.join(path.split('/')[:-1])
    data_name = path.split('/')[-1]
    print(path)
    print(data_path)
    print(data_name)
    ### data pathと保存pathを指定
    db = {
        'data_path': [data_path],
        'tiff_list': [data_name],
        'save_path0': './files',
        'save_folder': 'suite2p'
    }

    ops = {**default_ops(), **ops, **db}

    # save folderを指定
    save_folder = os.path.join(ops['save_path0'], ops['save_folder'])
    os.makedirs(save_folder, exist_ok=True)
    plane_folders = natsorted([
        f.path for f in os.scandir(save_folder) if f.is_dir() and f.name[:5]=='plane'])

    ops_path = [os.path.join(f, 'ops.npy') for f in plane_folders]

    if len(ops['h5py']):
        ops['input_format'] = 'h5'
    elif ops.get('mesoscan'):
        ops['input_format'] = 'mesoscan'
    elif not 'input_format' in ops:
        ops['input_format'] = 'tif'

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

    # plane_folders = natsorted([ f.path for f in os.scandir(save_folder) if f.is_dir() and f.name[:5]=='plane'])
    # ops_path = [os.path.join(f, 'ops.npy') for f in plane_folders][0]
    # ops = np.load(ops_path, allow_pickle=True).item()

    ops['images'] = ops['meanImg'].astype(np.uint8)

    return ops
