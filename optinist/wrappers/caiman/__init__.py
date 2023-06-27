from optinist.wrappers.caiman.cnmf import caiman_cnmf
from optinist.wrappers.caiman.motion_correction import caiman_mc

caiman_wrapper_dict = {
    "caiman": {
        "caiman_mc": {
            "function": caiman_mc,
            "conda_name": "caiman",
        },
        "caiman_cnmf": {
            "function": caiman_cnmf,
            "conda_name": "caiman",
        },
    }
}
