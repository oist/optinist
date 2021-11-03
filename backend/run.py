import os
import sys
sys.path.append('../optinist')
import json
import numpy as np
import pandas as pd
from PIL import Image
from collections import OrderedDict

from wrappers.data_wrapper import *


def run_code(wrapper_dict, flowList):
    print(flowList)
    print(wrapper_dict)

    info = {}
    inputs = None
    for item in flowList:
        print('-'*30)
        print('running')
        print(item)
        print('-'*30)
        if item.type == 'data':
            info[item.label] = {'path': ImageData(item.path, '')}
        elif item.type == 'algo':
            info[item.label] = wrapper_dict[item.label](*info[prev_label].values())

        prev_label = item.label

    results = OrderedDict()
    for item in flowList:
        results[item.label] = {}
        for k, v in info[item.label].items():
            results[item.label][k] = {}
            results[item.label][k]['path'] = v.path
            if type(v) is ImageData:
                print("ImageData") 
                results[item.label][k]['type'] = 'images'
                results[item.label][k]['max_index'] = len(v.data)
            elif type(v) is TimeSeriesData:
                print("TimeSeriesData")
                results[item.label][k]['type'] = 'timeseries'
            elif type(v) is CorrelationData:
                print("CorrelationData")
                results[item.label][k]['type'] = 'heatmap'

    print('results', results)

    return {'message': 'success', 'outputPaths': results}


if __name__ == '__main__':
    run_code()
