import sys
sys.path.append('../optinist')
from wrappers import wrapper_dict
from pytools.persistent_dict import PersistentDict

storage = PersistentDict("mystorage")

rule:
    input: 
        config["rules"]["dummy_image2scatter"]["input"]
    output: 
        touch(config["rules"]["dummy_image2scatter"]["output"])
    run:
        __func_config = config["rules"]["dummy_image2scatter"]
        info = storage.fetch(__func_config["input"])
        wrapper = dict2leaf(wrapper_dict, __func_config["path"].split('/'))
        info = wrapper["function"](*info.values())
        print(info)
        storage.store(__func_config["output"], info)


def dict2leaf(root_dict: dict, path_list):
    path = path_list.pop(0)
    if len(path_list) > 0:
        return dict2leaf(root_dict[path], path_list)
    else:
        return root_dict[path]