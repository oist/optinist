import os
import copy
import yaml

from wrappers.data_wrapper import *
from .utils import nest2dict, check_types
from .params import get_params


# def set_data(info, item, nwb_params):
#     if item['type'] == 'ImageFileNode':
#         filepath = os.path.join('..', 'optinist', 'config', f'nwb.yaml')
#         nwb_dict = copy.deepcopy(get_params(filepath))

#         if nwb_params != {}:
#             # 型を比較
#             nwb_dict = check_types(nest2dict(nwb_params), nwb_dict)

#         info = {'images': ImageData(item['data']['path'], '')}
#         nwb_dict['image_series']['external_file'] = info['images'].data
#         nwbfile = nwb_add_acquisition(nwb_dict)
#         nwbfile.create_processing_module(
#             name='ophys',
#             description='optical physiology processed data'
#         )
#         nwb_add_ophys(nwbfile)

#     elif item['type'] == 'CsvFileNode':
#         # import pdb; pdb.set_trace()
#         info = {item['key']: TimeSeriesData(item['data']['path'], '')}
#         nwbfile = None

#     info['nwbfile'] = nwbfile

#     return info


def set_data(dataType, nwb_params):
    if dataType == 'ImageFileNode':
        filepath = os.path.join('..', 'optinist', 'config', f'nwb.yaml')
        nwb_dict = copy.deepcopy(get_params(filepath))

        if nwb_params != {}:
            # 型を比較
            nwb_dict = check_types(nest2dict(nwb_params), nwb_dict)

        info = {'images': ImageData(item['data']['path'], '')}
        nwb_dict['image_series']['external_file'] = info['images'].data
        nwbfile = nwb_add_acquisition(nwb_dict)
        nwbfile.create_processing_module(
            name='ophys',
            description='optical physiology processed data'
        )
        nwb_add_ophys(nwbfile)

    elif dataType == 'CsvFileNode':
        info = {item['key']: TimeSeriesData(item['data']['path'], '')}
        nwbfile = None
    else:
        raise "error"

    info['nwbfile'] = nwbfile

    return info