from studio.app.optinist.wrappers.caiman.cnmf import caiman_cnmf
from studio.app.optinist.wrappers.caiman.cnmf_multisession import (
    caiman_cnmf_multisession,
)
from studio.app.optinist.wrappers.caiman.cnmfe import caiman_cnmfe
from studio.app.optinist.wrappers.caiman.motion_correction import caiman_mc

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
        "caiman_cnmfe": {
            "function": caiman_cnmfe,
            "conda_name": "caiman",
        },
        "cnmf_multisession": {
            "function": caiman_cnmf_multisession,
            "conda_name": "caiman",
        },
    }
}
