from .file_convert import suite2p_file_convert
from .registration import suite2p_registration
from .roi import suite2p_roi
from .spike_deconv import suite2p_spike_deconv


suite2p_wrapper_dict = {
    'suite2p': {
        'suite2p_file_convert': {
            'function': suite2p_file_convert,
        },
        'suite2p_registration': {
            'function': suite2p_registration,
        },
        'suite2p_roi': {
            'function': suite2p_roi,
        },
        'suite2p_spike_deconv': {
            'function': suite2p_spike_deconv
        },
    }
}
