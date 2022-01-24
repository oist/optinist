import os
import copy
import gc
from typing import List
from wrappers import wrapper_dict

from .params import get_params
from .utils import check_types


def run_algorithm(prev_info, item):

    if item['type'] == 'AlgorithmNode':

        filepath = os.path.join('..', 'optinist', 'config', f'{item["data"]["path"].split("/")[-1]}.yaml')
        params = copy.deepcopy(get_params(filepath))

        # parameterをint, floatに変換
        if 'param' in item['data'].keys() and item['data']['param'] is not None:
            params = check_types(item['data']['param'], params)

        wrapper = dict2leaf(
            wrapper_dict,
            item['data']['path'].split('/')
        )

        info = run_function(
            wrapper["function"],
            params,
            *prev_info.values(),
        )

        del wrapper, prev_info
        gc.collect()
    else:
        assert False, 'run_algorithm error'

    return info


def run_function(func_name, params, *args):
    info = func_name(params=params, *args)
    return info


def dict2leaf(root_dict: dict, path_list: List[str]):
    path = path_list.pop(0)
    if len(path_list) > 0:
        return dict2leaf(root_dict[path], path_list)
    else:
        return root_dict[path]
