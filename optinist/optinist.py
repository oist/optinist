

class Optinist:
	def __init__(self):
		print('init optinist')
		self._opts = None
	
	def set_data(self):
		print('load data')

	def set_params(self):
		print('set params')
		from algo import Caiman
		opts_dict = Caiman.set_params()
		print(opts_dict)
		return opts_dict

	def get_params(self):
		print('get params')

	def set_cluster(self):
		print('set cluster')

	def run_motion_correct(self):
		print('run motion correct')


if __name__ == '__main__':
	optinist = Optinist()
	optinist.set_params()
