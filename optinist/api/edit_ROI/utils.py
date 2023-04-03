import copy
import gc

from optinist.api.dir_path import DIRPATH
from optinist.api.edit_ROI import edit_roi_wrapper_dict


class EditRoiUtils:
    @classmethod
    def conda(cls, config):
        algo = config['algo']
        if "conda_yaml" in edit_roi_wrapper_dict[algo]:
            conda_yaml = edit_roi_wrapper_dict[algo]["conda_yaml"]
            return f"{DIRPATH.CONDAYML_DIR}/{conda_yaml}"

        return None

    @classmethod
    def excute(cls, config):
        algo = config['algo']
        node_dirpath = config['node_dirpath']
        action = config['action']
        params = config['params']

        func = copy.deepcopy(edit_roi_wrapper_dict[algo]["function"][action])
        output_info = func(node_dirpath, **params)

        del func
        gc.collect()
