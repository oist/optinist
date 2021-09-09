from wrappers.core import Wrapper

from caiman.source_extraction.cnmf.params import CNMFParams
from caiman.motion_correction import MotionCorrect
from caiman import load
from .set_cluster import set_cluster


class MotionCorrection(Wrapper):
	def __init__(self, env, opts=None):
		Wrapper.__init__(self, env)

		self.opts = CNMFParams()

		if opts is not None:
			self.opts.change_params(params_dict=opts)

	def get_params(self):
		return vars(self.opts)

	def run(self, info):
		info = self.env.run(info)
		mc = MotionCorrect(
			info['file_path'], dview=None, **self.opts.get_group('motion'))

		mc.motion_correct(save_movie=True)
		border_to_0 = 0 if mc.border_nan == 'copy' else mc.border_to_0

		info['mc'] = mc
		info['border_to_0'] = border_to_0

		return info

if __name__ == '__main__':
	import os
	env = MotionCorrection(None)
	file_path = os.path.join(
		'/Users', 'shogoakiyama', 'caiman_data', 
		'example_movies', 'Sue_2x_3000_40_-46.tif')
	env.run(file_path)
