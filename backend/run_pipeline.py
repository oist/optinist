import os
import sys
sys.path.append('../optinist')
import json
import numpy as np
import pandas as pd
from PIL import Image

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
            info[item.label] = wrapper_dict[item.label](*info[prev_label].values(), params=item.param)

        prev_label = item.label

    return info


if __name__ == '__main__':
    run_code()
