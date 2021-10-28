from wrappers.data_wrapper import *
from wrappers.args_check import args_check


@args_check
def suite2p_registration(image: ImageData, opts: dict=None):
    import numpy as np
    from suite2p import registration
    from suite2p import default_ops

    ######### REGISTRATION #########
    refImg = image.data
    if len(refImg.shape) == 3:
        refImg = refImg[0]

    if opts is None:
        ops = {**default_ops()}
    else:
        ops = opts

    ops = registration.register_binary(ops, refImg=refImg) # register binary

    # compute metrics for registration
    if ops.get('do_regmetrics', True) and ops['nframes']>=1500:
        ops = registration.get_pc_metrics(ops)

    ops['images'] = ops['refImg'].astype(np.uint8)
    ops['meanImgE'] = ops['meanImgE']

    return ops
