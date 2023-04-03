from .suite2p_edit_roi import (
    execute_add_ROI as suite2p_add,
    execute_merge_roi as suite2p_merge,
    excute_delete_roi as suite2p_delete,
)
from .lccd_edit_roi import (
    execute_add_ROI as lccd_add,
    execute_merge_roi as lccd_merge,
    excute_delete_roi as lccd_delete,
)


edit_roi_wrapper_dict = {
    'suite2p': {
        'conda_yaml': 'suite2p_env.yaml',
        'function': {
            'add': suite2p_add,
            'merge': suite2p_merge,
            'delete': suite2p_delete,
        },
    },
    'lccd': {
        'conda_yaml': 'lccd_env.yaml',
        'function': {
            'add': lccd_add,
            'merge': lccd_merge,
            'delete': lccd_delete,
        },
    },
}
