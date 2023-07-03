from studio.app.optinist.wrappers.suite2p.file_convert import suite2p_file_convert
from studio.app.optinist.wrappers.suite2p.registration import suite2p_registration
from studio.app.optinist.wrappers.suite2p.roi import suite2p_roi
from studio.app.optinist.wrappers.suite2p.spike_deconv import suite2p_spike_deconv

suite2p_wrapper_dict = {
    "suite2p": {
        "suite2p_file_convert": {
            "function": suite2p_file_convert,
            "conda_name": "suite2p",
        },
        "suite2p_registration": {
            "function": suite2p_registration,
            "conda_name": "suite2p",
        },
        "suite2p_roi": {
            "function": suite2p_roi,
            "conda_name": "suite2p",
        },
        "suite2p_spike_deconv": {
            "function": suite2p_spike_deconv,
            "conda_name": "suite2p",
        },
    }
}
