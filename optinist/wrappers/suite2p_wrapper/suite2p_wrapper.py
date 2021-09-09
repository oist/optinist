import os
import yaml

from base import Base

class Suite2p(Base):
	def __init__(self):
		from suite2p import default_ops

		self.fnames = None
		self.opts = default_ops()

	def set_data(self, file_path):
		self.fnames = [file_path]

	def set_params(self, opts_dict=None):

		if opts_dict is None:

			with open(os.path.join('./algo', 'config', 'suite2p.yaml')) as f:
				opts_dict = yaml.load(f, Loader=yaml.FullLoader)
				assert opts_dict

		for k, v in opts_dict.items():
			self.opts[k] = v

	def get_params(self):
		return self.opts

	def run(self):
		import suite2p
		data_path = []
		tiff_list = []
		for fname in self.fnames:
			data_path.append('/'.join(fname.split('/')[:-1]))
			tiff_list.append(fname.split('/')[-1])

		db = {
			'data_path': data_path,
			'save_dir': './logs',
			'save_path0': './logs',
			'tiff_list': tiff_list
		}

		output_ops = suite2p.run_s2p(ops=self.opts, db=db)
		print('output path: ', output_ops['save_path'])

	def run_motion_correct(self):
		from suite2p.registration import register
		refImg = register.pick_initial_reference(ops)

	def plot_output(self):
		import os
		import numpy as np
		import matplotlib.pyplot as plt

		path = os.path.join('./logs', 'suite2p', 'plane0', 'ops.npy')

		assert os.path.exists(path)

		output_op = np.load(path, allow_pickle=True).item()

		plt.figure(figsize=(16, 4))
		plt.subplot(1, 4, 1)
		plt.imshow(output_op['refImg'], cmap='gray')
		plt.title("Reference Image for Registration")

		plt.subplot(1, 4, 2)
		plt.imshow(output_op['max_proj'], cmap='gray')
		plt.title("Registered Image, Max Projection")

		plt.subplot(1, 4, 3)
		plt.imshow(output_op['meanImg'], cmap='gray')
		plt.title("Mean registered image")

		plt.subplot(1, 4, 4)
		plt.imshow(output_op['meanImgE'], cmap='gray')
		plt.title("High-pass filtered Mean registered image")

		plt.show()

if __name__ == '__main__':
	algo = Suite2p()
	file_path = os.path.join(
		'./data', 'Sue_2x_3000_40_-46.tif')
	algo.set_data(file_path)
	# print(algo.get_params())
	algo.set_params()
	# print(algo.get_params())
	# algo.set_cluster()
	algo.run()
	algo.plot_output()
