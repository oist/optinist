from studio.core.edit_ROI.wrappers.caiman_edit_roi import (
    excute_delete_roi as caiman_delete,
)
from studio.core.edit_ROI.wrappers.caiman_edit_roi import execute_add_ROI as caiman_add
from studio.core.edit_ROI.wrappers.caiman_edit_roi import (
    execute_merge_roi as caiman_merge,
)
from studio.core.edit_ROI.wrappers.lccd_edit_roi import excute_delete_roi as lccd_delete
from studio.core.edit_ROI.wrappers.lccd_edit_roi import execute_add_ROI as lccd_add
from studio.core.edit_ROI.wrappers.lccd_edit_roi import execute_merge_roi as lccd_merge
from studio.core.edit_ROI.wrappers.suite2p_edit_roi import (
    excute_delete_roi as suite2p_delete,
)
from studio.core.edit_ROI.wrappers.suite2p_edit_roi import (
    execute_add_ROI as suite2p_add,
)
from studio.core.edit_ROI.wrappers.suite2p_edit_roi import (
    execute_merge_roi as suite2p_merge,
)

edit_roi_wrapper_dict = {
    "suite2p": {
        "conda_name": "suite2p",
        "function": {
            "add": suite2p_add,
            "merge": suite2p_merge,
            "delete": suite2p_delete,
        },
    },
    "lccd": {
        "conda_name": None,
        "function": {
            "add": lccd_add,
            "merge": lccd_merge,
            "delete": lccd_delete,
        },
    },
    "caiman": {
        "conda_name": None,
        "function": {
            "add": caiman_add,
            "merge": caiman_merge,
            "delete": caiman_delete,
        },
    },
}
