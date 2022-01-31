from wrappers import wrapper_dict
from pytools.persistent_dict import PersistentDict
import traceback

storage = PersistentDict("mystorage")

rule:
    input:
        config["rules"]["dummy_image8image2image8time"]["input"]
    output:
        touch(config["rules"]["dummy_image8image2image8time"]["output"])
    run:
        __func_config = config["rules"]["dummy_image8image2image8time"]
        try:
            input_files = __func_config["input"]
            info = {}
            for path in input_files:
                info.update(storage.fetch(path))

            params = __func_config["params"]
            wrapper = dict2leaf(wrapper_dict, __func_config["path"].split('/'))

            output_info = wrapper["function"](params=params, *info.values())
            storage.store(__func_config["output"], output_info)
        except Exception as e:
            error_message  = list(traceback.TracebackException.from_exception(e).format())[-1]
            storage.store(__func_config["output"], error_message)


def dict2leaf(root_dict: dict, path_list):
    path = path_list.pop(0)
    if len(path_list) > 0:
        return dict2leaf(root_dict[path], path_list)
    else:
        return root_dict[path]
