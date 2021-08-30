class Base:
	def __init__(self):
		print('init optinist')
		self._opts = None
	
	def set_data(self):
		print('load data')

	def set_params(self):
		print('set params')

	def set_cluster(self):
		print('set cluster')

	def run_motion_correct(self):
		print('run motion correct')

	def get_params(self):
		print('get params')
