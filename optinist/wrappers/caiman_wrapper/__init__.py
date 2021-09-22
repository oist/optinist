from .motion_correction import caiman_mc
from .cnmf import caiman_cnmf
from .plot_contours import plot_contours_nb

caiman_wrapper_dict = {
	'caiman_mc': caiman_mc,
	'caiman_cnmf': caiman_cnmf,
	'plot_contours_nb': plot_contours_nb
}