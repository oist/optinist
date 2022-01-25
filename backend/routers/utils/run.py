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

        prev_info = order_args(wrapper, prev_info)

        info = run_function(
            wrapper["function"],
            params,
            *prev_info.values(),
            # prev_info,
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

def order_args(wrapper, prev_info):
    import inspect
    sig = inspect.signature(wrapper["function"])

    # 引数名の順番に揃える
    new_args = {}
    for key in sig.parameters.keys():
        if key in prev_info.keys():
            new_args[key] = prev_info[key]

    for key in prev_info.keys():
        if not key in new_args:
            new_args[key] = prev_info[key]

    return new_args
