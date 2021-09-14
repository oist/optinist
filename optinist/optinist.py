import os
import wrappers.caiman_wrapper as cm
import wrappers.suite2p_wrapper as s2p


class Optinist(object):
	def __init__(self, opt=None):
		self.opt = opt

	def run(self, info):
		return info


def create_env():
	env = Optinist()
	env = cm.MotionCorrection(env)
	env = cm.MemoryMapping(env)
	env = cm.CNMF(env)
	return env


if __name__ == '__main__':

	file_path = os.path.join(
		'/Users', 'shogoakiyama', 'caiman_data', 
		'example_movies', 'Sue_2x_3000_40_-46.tif')

	env = create_env()
	info = {'file_path': file_path}
	env.run(info)
