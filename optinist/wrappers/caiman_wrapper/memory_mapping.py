from wrappers.core import Wrapper

from caiman.source_extraction.cnmf.params import CNMFParams
from caiman import save_memmap, load_memmap

import numpy as np


class MemoryMapping(Wrapper):
	def __init__(self, env, opts=None):
		Wrapper.__init__(self, env)
		print(self.env)

	def run(self, info):
		info = self.env.run(info)

		fname_new = save_memmap(
			info['mc'].mmap_file, base_name='memmap_', order='C',border_to_0=info['border_to_0'])

		# now load the file
		Yr, dims, T = load_memmap(fname_new)
		images = np.reshape(Yr.T, [T] + list(dims), order='F') 

		info['images'] = images

		return info

if __name__ == '__main__':
	import os
	env = MotionCorrection(None)
	file_path = os.path.join(
		'/Users', 'shogoakiyama', 'caiman_data', 
		'example_movies', 'Sue_2x_3000_40_-46.tif')
	env.run(file_path)
