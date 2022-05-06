from fastapi import APIRouter
from typing import List
import inspect

from optinist.routers.const import NOT_DISPLAY_ARGS_LIST
from optinist.routers.model import Algo, Arg, Return
from optinist.wrappers import wrapper_dict

router = APIRouter()


def create_args_list(args_dict):
    return [
        Arg(
            name=x.name,
            type=x.annotation.__name__,
            isNone=x.default is None,
        )
        for x in args_dict
        if x.name not in NOT_DISPLAY_ARGS_LIST
    ]


def create_return_list(return_dict):
    return [
        Return(
            name=k,
            type=v.__name__
        )
        for k, v in return_dict
    ]


def create_parent_key(parent_key, key):
    if parent_key == '':
        return key
    else:
        return f'{parent_key}/{key}'


def get_nest_dict(parent_value, parent_key):
    algo_dict = {}
    for key, value in parent_value.items():
        algo_dict[key] = {}
        if isinstance(value, dict) and 'function' not in value:
            algo_dict[key]['children'] = get_nest_dict(
                value,
                create_parent_key(parent_key, key)
            )
        else:
            sig = inspect.signature(value['function'])
            returns_list = None
            if sig.return_annotation is not inspect._empty:
                returns_list = create_return_list(sig.return_annotation.items())

            algo_dict[key] = Algo(
                args=create_args_list(sig.parameters.values()),
                returns=returns_list,
                parameter=value['parameter'] if 'parameter' in value else None,
                path=create_parent_key(parent_key, key),
            )

    return algo_dict


@router.get("/algolist")
async def run() -> List:
    {
        'caiman': {
            'children': {
                'caiman_mc' : {
                    'args': ['images', 'timeseries'],
                    'return': ['images'],
                    'path': 'caiman/caiman_mc'
                },
                'caiman_cnmf': {
                    'args': ['images', 'timeseries'],
                    'return': ['images'],
                    'path': 'caiman/caiman_mc'
                }
            }
        }
    }
    return get_nest_dict(wrapper_dict, '')
