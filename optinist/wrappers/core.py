class Wrapper(object):

	def __init__(self, env):
		self.env = env

	@classmethod
	def run(self):
		return self.env.run()
