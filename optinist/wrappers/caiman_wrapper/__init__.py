from .motion_correction import caiman_mc
from .cnmf import caiman_cnmf


caiman_wrapper_dict = {
    'caiman': {
        'caiman_mc': { 
            'function': caiman_mc,
            'parameter': 'caiman_mc.yaml'
        },
        'caiman_cnmf': {
            'function': caiman_cnmf,
            'parameter': 'caiman_cnmf.yaml'
        },
    }
}
