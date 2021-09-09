from wrappers.core import Wrapper
from suite2p import default_ops
from suite2p.registration import register


class MotionCorrection(Wrapper):
	def __init__(self, env, opts=None):
		Wrapper.__init__(self, env)
		self.ops = default_ops()

	def get_params(self):
		return vars(self.ops)

	def run(self, info):
		info = self.env.run(info)
		refImg = register.pick_initial_reference(self.ops)

		return info

if __name__ == '__main__':
	import os
	env = MotionCorrection(None)
	file_path = os.path.join(
		'/Users', 'shogoakiyama', 'caiman_data', 
		'example_movies', 'Sue_2x_3000_40_-46.tif')
	env.run(file_path)
