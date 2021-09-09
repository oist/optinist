from wrappers.core import Wrapper

from caiman.source_extraction.cnmf.params import CNMFParams
from caiman.source_extraction.cnmf import cnmf

import numpy as np


class CNMF(Wrapper):
	def __init__(self, env, opts=None):
		Wrapper.__init__(self, env)
		print(self.env)

	def run(self, info):
		info = self.env.run(info)

		cnm = cnmf.CNMF(n_processes=0)
		cnm = cnm.fit(info['images'])

		return info

if __name__ == '__main__':
	import os
	env = MotionCorrection(None)
	file_path = os.path.join(
		'/Users', 'shogoakiyama', 'caiman_data', 
		'example_movies', 'Sue_2x_3000_40_-46.tif')
	env.run(file_path)
