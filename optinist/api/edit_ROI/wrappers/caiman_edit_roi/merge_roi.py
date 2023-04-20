from optinist.api.dataclass.dataclass import *


def execute_merge_roi(node_dirpath: str, ids: list):
    import numpy as np

    from .utils import get_roi

    # load data
    cnmf_data = np.load(f'{node_dirpath}/caiman_cnmf.npy', allow_pickle=True).item()
    is_cell = cnmf_data.get('is_cell')
    ims = cnmf_data.get('ims')
    fluorescence = cnmf_data.get('fluorescence')

    # get merging ROI
    merging_ROIs = []
    [merging_ROIs.append(ims[id, :, :]) for id in ids]

    # get merged ROI
    merged_ROI = np.maximum.reduce(merging_ROIs)
    is_cell[ids] = False
    is_cell = np.append(is_cell, True)
    ims = np.concatenate((ims, merged_ROI[np.newaxis, :, :]), axis=0)

    # get merged F
    merged_f = np.mean((fluorescence[ids, :]), axis=0)
    fluorescence = np.vstack([fluorescence, merged_f])

    cell_roi = get_roi(ims)

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
