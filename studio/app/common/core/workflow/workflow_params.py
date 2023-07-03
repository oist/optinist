from studio.app.common.core.utils.config_handler import ConfigReader
from studio.app.common.core.utils.filepath_finder import find_param_filepath


def get_typecheck_params(message_params, name):
    default_params = ConfigReader.read(find_param_filepath(name))
    if message_params != {} and message_params is not None:
        return check_types(nest2dict(message_params), default_params)
    return default_params


def check_types(params, default_params):
    for key in params.keys():
        if isinstance(params[key], dict):
            params[key] = check_types(params[key], default_params[key])
        else:
            if not isinstance(type(params[key]), type(default_params[key])):
                data_type = type(default_params[key])
                p = params[key]
                if isinstance(data_type, str):
                    params[key] = str(p)
                elif isinstance(data_type, float):
                    params[key] = float(p)
                elif isinstance(data_type, int):
                    params[key] = int(p)

    return params


def nest2dict(value):
    nwb_dict = {}
    for _k, _v in value.items():
        if _v["type"] == "child":
            nwb_dict[_k] = _v["value"]
        elif _v["type"] == "parent":
            nwb_dict[_k] = nest2dict(_v["children"])

    return nwb_dict
