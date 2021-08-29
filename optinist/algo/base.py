class Base:
	def __init__(self):
		print('init optinist')
		self._params = None

	def run(self) -> list:
		print('run')

	def set_data(self, data_path: str):
		print('load data')

	def set_params(self, params: dict):
		print('set params')
		self.params = params

	def get_params(self) -> dict:
		print('get params')
		return self.params

	def set_cluster(self):
		print('set cluster')

	def run_motion_correct(self):
		print('run motion correct')


if __name__ == '__main__':
	algo = Base()
	algo.set_params({'a': 1, 'b': 2})
	print(algo.get_params())
	# algo.set_params(['a'])
	print(algo.get_params())
	algo.run()
