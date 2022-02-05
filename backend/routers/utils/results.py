from collections import OrderedDict
from wrappers.data_wrapper import *


def get_results(info, path):
    results = OrderedDict()
    results[path] = {}

    for k, v in info.items():
        if isinstance(v, ImageData):
            results[path][k] = {}
            results[path][k]['path'] = v.json_path
            results[path][k]['type'] = 'images'
            if v.data.ndim == 3:
                results[path][k]['max_index'] = len(v.data)
            elif v.data.ndim == 2:
                results[path][k]['max_index'] = 1
        elif isinstance(v, TimeSeriesData):
            results[path][k] = {}
            results[path][k]['path'] = v.path
            results[path][k]['type'] = 'timeseries'
            results[path][k]['max_index'] = len(v.data)
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
        elif isinstance(v, BarData):
            results[path][k] = {}
            results[path][k]['path'] = v.path
            results[path][k]['type'] = 'bar'
        else:
            pass

    return results