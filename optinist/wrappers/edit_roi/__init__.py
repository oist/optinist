import os
import gc
import copy

from .suite2p_edit_roi import (
    execute_add_ROI as suite2p_add,
    execute_merge_roi as suite2p_merge,
    excute_delete_roi as suite2p_delete,
)
from fastapi import HTTPException, status
from optinist.api.dir_path import DIRPATH


edit_roi_wrapper_dict = {
    'suite2p': {
        'conda_yaml': 'suite2p_env.yaml',
        'function': {
            'add': suite2p_add,
            'merge': suite2p_merge,
            'delete': suite2p_delete,
        },
    },
    'caiman': {
        'conda_yaml': 'caiman_env.yaml',
        'function': {
            # 'add': execute_add_ROI,
            # 'merge': execute_merge_roi,
            # 'delete': excute_delete_roi
        },
    },
}


class EditRoiUtils:
    @classmethod
    def params(cls, config):
        return config['params']

    @classmethod
    def conda(cls, config):
        algo = cls.algo(config)
        if not algo:
            return None
        print(algo)
        if "conda_yaml" in edit_roi_wrapper_dict[algo]:
            conda_yaml = edit_roi_wrapper_dict[algo]["conda_yaml"]
            return f"{DIRPATH.CONDAYML_DIR}/{conda_yaml}"

        return None

    @classmethod
    def algo(cls, config):
        if 'filepath' not in config:
            return None
        filepath = config['filepath']
        if (
            not os.path.exists(filepath)
            and os.path.commonpath([DIRPATH.OPTINIST_DIR, filepath])
            != DIRPATH.OPTINIST_DIR
        ):
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND)
        if os.path.basename(filepath) != 'cell_roi.json':
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST)

        algo_list = edit_roi_wrapper_dict.keys()
        return next((algo for algo in algo_list if algo in filepath), None)

    @classmethod
    def excute(cls, config):
        # import suite2p
        algo = cls.algo(config)
        action = config['action']

        node_dirpath = os.path.dirname(config['filepath'])
        params = cls.params(config)
        print(algo, node_dirpath, params)
        func = copy.deepcopy(edit_roi_wrapper_dict[algo]["function"][action])
        print(func.__name__)
        output_info = func(node_dirpath, **params)

        del func
        gc.collect()
