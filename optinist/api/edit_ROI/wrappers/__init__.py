from .suite2p_edit_roi import (
    execute_add_ROI as suite2p_add,
    execute_merge_roi as suite2p_merge,
    excute_delete_roi as suite2p_delete,
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
}
