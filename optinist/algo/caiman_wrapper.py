import os
import caiman as cm

from base import Base

class Caiman(Base):
	def __init__(self):
		from caiman.source_extraction.cnmf import params as params

		self.fnames = None
		self.dview = None
		self.opts = params.CNMFParams()
	
	def set_data(self, file_path):
		self.fnames = [file_path]

	# @classmethod
	def set_params(self, opts_dict=None):
		
		if opts_dict is None:
			import yaml
			import os

			with open(os.path.join('./algo', 'config', 'caiman.yaml')) as f:
					opts_dict = yaml.load(f, Loader=yaml.FullLoader)
					# print(opts_dict)
					assert opts_dict

		self.opts.change_params(params_dict=opts_dict)
		# import pdb; pdb.set_trace()

	def get_params(self):
		return vars(self.opts)

	def set_cluster(self):
		print('set cluster')
		import caiman as cm
		import platform

		# set cluster
		if 'dview' in locals():
			print('stop cluster')
			cm.stop_server(dview=self.dview)

		single_thread = False
		if platform.system() == 'Darwin':
			# bug issue in [https://github.com/flatironinstitute/CaImAn/issues/206]
			# and check code [https://github.com/flatironinstitute/CaImAn/blob/master/caiman/cluster.py#L414]
			single_thread = True

		c, dview, n_processes = cm.cluster.setup_cluster(
			backend='local', n_processes=None, single_thread=single_thread)

		self.dview = dview

	def run_motion_correct(self):
		print('run motion correct')
		assert self.fnames
		assert self.opts

		from caiman.motion_correction import MotionCorrect

		mc = MotionCorrect(self.fnames, dview=self.dview, **self.opts.get_group('motion'))

		mc.motion_correct(save_movie=True)
		m_els = cm.load(mc.fname_tot_els)
		border_to_0 = 0 if mc.border_nan == 'copy' else mc.border_to_0

if __name__ == '__main__':
	algo = Caiman()
	file_path = os.path.join(
		'/Users', 'shogoakiyama', 'caiman_data', 
		'example_movies', 'Sue_2x_3000_40_-46.tif')
	algo.set_data(file_path)
	# print(algo.get_params())
	algo.set_params()
	# print(algo.get_params())
	algo.set_cluster()
	algo.run_motion_correct()
