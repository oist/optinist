from optinist.api.edit_ROI.wrappers.caiman_edit_roi import (
    excute_delete_roi as caiman_delete,
)
from optinist.api.edit_ROI.wrappers.caiman_edit_roi import execute_add_ROI as caiman_add
from optinist.api.edit_ROI.wrappers.caiman_edit_roi import (
    execute_merge_roi as caiman_merge,
)
from optinist.api.edit_ROI.wrappers.lccd_edit_roi import (
    excute_delete_roi as lccd_delete,
)
from optinist.api.edit_ROI.wrappers.lccd_edit_roi import execute_add_ROI as lccd_add
from optinist.api.edit_ROI.wrappers.lccd_edit_roi import execute_merge_roi as lccd_merge
from optinist.api.edit_ROI.wrappers.suite2p_edit_roi import (
    excute_delete_roi as suite2p_delete,
)
from optinist.api.edit_ROI.wrappers.suite2p_edit_roi import (
    execute_add_ROI as suite2p_add,
)
from optinist.api.edit_ROI.wrappers.suite2p_edit_roi import (
    execute_merge_roi as suite2p_merge,
)

edit_roi_wrapper_dict = {
    "suite2p": {
        "conda_yaml": "suite2p_env.yaml",
        "function": {
            "add": suite2p_add,
            "merge": suite2p_merge,
            "delete": suite2p_delete,
        },
    },
    "lccd": {
        "conda_yaml": None,
        "function": {
            "add": lccd_add,
            "merge": lccd_merge,
            "delete": lccd_delete,
        },
    },
    "caiman": {
        "conda_yaml": None,
        "function": {
            "add": caiman_add,
            "merge": caiman_merge,
            "delete": caiman_delete,
        },
    },
}
