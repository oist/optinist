import numpy as np
import caiman


def plot_contours_nb(images, cnm):
	Cn = caiman.local_correlations(images.transpose(1, 2, 0))
	Cn[np.isnan(Cn)] = 0
	cnm.estimates.plot_contours_nb(img=Cn)


if __name__ == '__main__':
	import os
	from motion_correction import caiman_mc
	from cnmf import caiman_cnmf
	import caiman

	info = {}
	file_path = os.path.join(
		'/Users', 'shogoakiyama', 'caiman_data', 
		'example_movies', 'Sue_2x_3000_40_-46.tif')
	info['caiman_mc'] = caiman_mc(file_path)
	info['caiman_cnmf'] = caiman_cnmf(info['caiman_mc']['images'])
	plot_contours_nb(info['caiman_mc']['images'], info['caiman_cnmf']['cnm'])
