from wrappers import wrapper_dict
from pytools.persistent_dict import PersistentDict
import traceback
import os
import pickle

storage = PersistentDict("mystorage")


def dict2leaf(root_dict: dict, path_list):
    path = path_list.pop(0)
    if len(path_list) > 0:
        return dict2leaf(root_dict[path], path_list)
    else:
        return root_dict[path]


def change_dict_key_exist(d, old_key, new_key):
    if old_key in d:
        d[new_key] = d.pop(old_key)


def run_script(__func_config):
    try:
        input_files = __func_config["input"]
        return_arg = __func_config["return_arg"]
        info = {}

        for path in input_files:
            # info.update(storage.fetch(path))
            # with open(path, 'rb') as f:
            #     info.update(pickle.load(f))
            path = path.split(".")[0] + ".pkl"
            with open(path, 'rb') as f:
                data = pickle.load(f)
                info.update(data)

        params = __func_config["params"]
        wrapper = dict2leaf(wrapper_dict, __func_config["path"].split('/'))
        print(wrapper)

        for return_name, arg_name in return_arg.items():
            change_dict_key_exist(info, return_name, arg_name)

        for key in list(info):
            if key != "nwbfile" and key not in return_arg.values():
                info.pop(key)

        output_info = wrapper["function"](params=params, **info)
        storage.store(__func_config["output"], output_info)

        outdir = ".".join(__func_config["output"].split("/")[-1:])
        # os.makedirs(outdir, exist_ok=True)
        with open(__func_config["output"], 'wb') as f:
            pickle.dump(output_info, f)
        
        print("output: ", __func_config["output"])
        
    except Exception as e:
        error_message  = list(traceback.TracebackException.from_exception(e).format())[-2:]
        # storage.store(__func_config["output"], error_message)
        with open(__func_config["output"], 'wb') as f:
            pickle.dump(error_message, f)
        raise "error"
