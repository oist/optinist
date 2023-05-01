from optinist.wrappers.suite2p_wrapper.file_convert import suite2p_file_convert
from optinist.wrappers.suite2p_wrapper.registration import suite2p_registration
from optinist.wrappers.suite2p_wrapper.roi import suite2p_roi
from optinist.wrappers.suite2p_wrapper.spike_deconv import suite2p_spike_deconv

suite2p_wrapper_dict = {
    "suite2p": {
        "suite2p_file_convert": {
            "function": suite2p_file_convert,
            "conda_name": "suite2p",
            "conda_yaml": "suite2p_env.yaml",
        },
        "suite2p_registration": {
            "function": suite2p_registration,
            "conda_name": "suite2p",
            "conda_yaml": "suite2p_env.yaml",
        },
        "suite2p_roi": {
            "function": suite2p_roi,
            "conda_name": "suite2p",
            "conda_yaml": "suite2p_env.yaml",
        },
        "suite2p_spike_deconv": {
            "function": suite2p_spike_deconv,
            "conda_name": "suite2p",
            "conda_yaml": "suite2p_env.yaml",
        },
    }
}
