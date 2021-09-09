from caiman import stop_server
from caiman.cluster import setup_cluster
import platform


def set_cluster():
	print('set cluster')

	single_thread = False
	if platform.system() == 'Darwin':
		# bug issue in [https://github.com/flatironinstitute/CaImAn/issues/206]
		# and check code [https://github.com/flatironinstitute/CaImAn/blob/master/caiman/cluster.py#L414]
		single_thread = True

	c, dview, n_processes = setup_cluster(
		backend='local', n_processes=None, single_thread=single_thread)