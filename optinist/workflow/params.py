import os
import copy
import yaml
from optinist.cui_api.config_reader import ConfigReader


def get_typecheck_params(message_params, name):
    params = ConfigReader.read(name)
    if message_params != {} and message_params is not None:
        params = check_types(nest2dict(message_params), params)
    return params


def check_types(params, default_params):
    for key in params.keys():
        if isinstance(params[key], dict):
            params[key] = check_types(params[key], default_params[key])
        else:
            if type(params[key]) != type(default_params[key]):
                data_type = type(default_params[key])
                p = params[key]
                if data_type == str:
                    params[key] = str(p)
                elif data_type == float:
                    params[key] = float(p)
                elif data_type == int:
                    params[key] = int(p)

    return params


def nest2dict(value):
    nwb_dict = {}
    for _k, _v in value.items():
        if _v['type'] == 'child':
            nwb_dict[_k] = _v['value']
        elif _v['type'] == 'parent':
            nwb_dict[_k] = nest2dict(_v['children'])

    return nwb_dict
