


def caiman_cnmf(images):
	from caiman.source_extraction.cnmf import cnmf
	info = {}

	cnm = cnmf.CNMF(n_processes=0)
	cnm = cnm.fit(images)
	
	info['cnm'] = cnm

	return info

if __name__ == '__main__':
	import os
	import numpy as np
	from motion_correction import caiman_mc
	import caiman

	info = {}
	file_path = os.path.join(
		'/Users', 'shogoakiyama', 'caiman_data', 
		'example_movies', 'Sue_2x_3000_40_-46.tif')
	info['caiman_mc'] = caiman_mc(file_path)
	info['caiman_cnmf'] = caiman_cnmf(info['caiman_mc']['images'])
