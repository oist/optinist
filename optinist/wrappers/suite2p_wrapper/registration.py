from wrappers.data_wrapper import *
from wrappers.args_check import args_check


@args_check
def suite2p_registration(image: ImageData, ops: Suite2pData=None):
    refImg = image.data
    ops = ops.data

    import numpy as np
    from suite2p import registration
    from suite2p import default_ops

    ######### REGISTRATION #########
    if len(refImg.shape) == 3:
        refImg = refImg[0]

    if ops is None:
        ops = {**default_ops()}

    ops = registration.register_binary(ops, refImg=refImg) # register binary

    # compute metrics for registration
    if ops.get('do_regmetrics', True) and ops['nframes']>=1500:
        ops = registration.get_pc_metrics(ops)

    info = {}
    info['images'] = ImageData(ops['refImg'].astype(np.uint8), 'refImg')
    info['meanImgE'] = ImageData(ops['meanImgE'], 'meanImgE')
    info['ops'] = Suite2pData(ops)

    return info
