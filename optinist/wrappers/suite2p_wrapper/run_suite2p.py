from wrappers.data_wrapper import *
from wrappers.args_check import args_check


@args_check
def run_suite2p(image: ImageData, ops: dict=None)
    file_path = image.path
    data_path = '/'.join(file_path.split('/')[:-1])
    data_name = file_path.split('/')[-1]

    db = {
        'data_path': [data_path],
        'tiff_list': [data_name],
        'save_path0': './files',
        'save_folder': 'suite2p'
    }

    if opts is None:
        ops = {**default_ops(), **db}
    else:
        ops = {**ops, **db}

	output_ops = suite2p.run_s2p(ops=ops, db=db)

	return output_ops
