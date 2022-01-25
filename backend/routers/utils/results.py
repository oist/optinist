from collections import OrderedDict
from wrappers.data_wrapper import *


def get_results(info, item):
    results = OrderedDict()
    path = item['data']['path']
    results[path] = {}

    for k, v in info.items():
        if isinstance(v, ImageData):
            results[path][k] = {}
            results[path][k]['path'] = v.json_path
            results[path][k]['type'] = 'images'
            results[path][k]['max_index'] = len(v.data)
        elif isinstance(v, TimeSeriesData):
            results[path][k] = {}
            results[path][k]['path'] = v.path
            results[path][k]['type'] = 'timeseries'
        elif isinstance(v, CorrelationData):
            results[path][k] = {}
            results[path][k]['path'] = v.path
            results[path][k]['type'] = 'heatmap'
        elif isinstance(v, RoiData):
            results[path][k] = {}
            results[path][k]['path'] = v.path
            results[path][k]['type'] = 'roi'
        elif isinstance(v, ScatterData):
            results[path][k] = {}
            results[path][k]['path'] = v.path
            results[path][k]['type'] = 'scatter'
        else:
            pass

    return results