from optinist.api.dataclass.dataclass import *


def execute_add_ROI(node_dirpath, posx, posy, sizex, sizey):
    import numpy as np

    from .utils import create_mask, get_roi

    # load data
    cnmf_data = np.load(f'{node_dirpath}/caiman_cnmf.npy', allow_pickle=True).item()
    dims = cnmf_data['dims']
    is_cell = cnmf_data.get('is_cell')
    ims = cnmf_data.get('ims')
    images = cnmf_data.get('images')
    fluorescence = cnmf_data.get('fluorescence')

    # Create mask for new roi
    new_roi = create_mask(posx, posy, sizex, sizey, dims)
    ims = np.concatenate((ims, new_roi[np.newaxis, :, :]), axis=0)
    is_cell = np.append(is_cell, True)

    # extract fluorescence of new ROI
    num_frames = images.shape[0]
    reshapeImages = images.reshape([images.shape[2] * images.shape[1], num_frames])
    new_fluorescence = np.zeros(num_frames)
    new_fluorescence = np.mean(reshapeImages[new_roi.reshape(-1, 1)[:, 0] > 0], axis=0)
    fluorescence = np.vstack([fluorescence, new_fluorescence])

    cell_roi = get_roi(ims)

    # save data
    cnmf_data['ims'] = ims
    cnmf_data['is_cell'] = is_cell
    cnmf_data['fluorescence'] = fluorescence

    info = {
        'fluorescence': FluoData(fluorescence, file_name='fluorescence'),
        'cell_roi': RoiData(np.nanmax(cell_roi[is_cell], axis=0), file_name='cell_roi'),
        'cnmf_data': CaimanCnmfData(cnmf_data),
    }

    for v in info.values():
        if isinstance(v, BaseData):
            v.save_json(node_dirpath)
