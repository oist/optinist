from .motion_correction import caiman_mc
from .cnmf import caiman_cnmf


caiman_wrapper_dict = {
	'caiman': {
		'caiman_mc': caiman_mc,
		'caiman_cnmf': caiman_cnmf,
	}
}
