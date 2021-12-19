from fastapi import APIRouter
from typing import List
import inspect

from wrappers import wrapper_dict

router = APIRouter()

def get_nest_dict(value, parent_k):
    algo_dict = {}
    for _k, _v in value.items():
        algo_dict[_k] = {}
        if type(_v) is dict and 'function' not in _v.keys():
            algo_dict[_k]['children'] = get_nest_dict(
                _v, parent_k+'/'+_k if parent_k != '' else _k)
        else:
            # get args
            sig = inspect.signature(_v['function'])
            algo_dict[_k]['args'] = [
                {
                    'name': x.name, 
                    'type': x.annotation.__name__
                }
                for x in sig.parameters.values()
                if x.name != 'params'
            ]

            # get returns
            if sig.return_annotation is not inspect._empty:
                algo_dict[_k]['returns'] = [
                    {
                        'name': k,
                        'type': v.__name__
                    }
                    for k, v in sig.return_annotation.items()
                ]
            
            # parameter path
            if 'parameter' in _v.keys():
                algo_dict[_k]['parameter'] = _v['parameter']
            else:
                algo_dict[_k]['parameter'] = ''
            
            # path
            algo_dict[_k]['path'] = parent_k + '/' + _k

    return algo_dict


@router.get("/algolist")
async def run() -> List:
    # print(wrapper_dict.keys())
    {
        'caiman': {
            'children': {
                'caiman_mc' : {
                    'args': ['images', 'timeseries'],
                    'path': 'caiman/caiman_mc'
                },
                'caiman_cnmf': {
                    'args': ['images', 'timeseries'],
                    'path': 'caiman/caiman_mc'
                }
            }
        }
    }
    {
        'caiman.caiman_mc': {

        },
        'caiman.caiman_cnmf': {

        },
    }

    algo_dict = get_nest_dict(wrapper_dict, '')

    return algo_dict
