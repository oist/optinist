from wrappers.data_wrapper import *
from wrappers.args_check import args_check


@args_check
def suite2p_registration(
        ops: Suite2pData=None, params: dict=None
    ) -> {'ops': Suite2pData, 'images': ImageData}:
    # refImg = image.data
    ops = ops.data
    refImg = ops['meanImg']

    import numpy as np
    from suite2p import registration, io, default_ops

    ######### REGISTRATION #########
    if len(refImg.shape) == 3:
        refImg = refImg[0]

    # if ops is None:
    #     ops = {**default_ops()}

    ops = {**ops, **params}

    ops = registration.register_binary(ops, refImg=refImg) # register binary

    # compute metrics for registration
    if ops.get('do_regmetrics', True) and ops['nframes']>=1500:
        ops = registration.get_pc_metrics(ops)

    info = {}
    from suite2p import io
    images = io.BinaryFile(read_filename=ops['reg_file'], Ly=ops['Ly'], Lx=ops['Lx']).data
    info['images'] = ImageData(images, func_name='suite2p_registration', file_name='images')
    info['refImg'] = ImageData(ops['refImg'].astype(np.uint8), func_name='suite2p_registration', file_name='refImg')
    info['meanImgE'] = ImageData(ops['meanImgE'].astype(np.uint8), func_name='suite2p_registration', file_name='meanImgE')
    info['ops'] = Suite2pData(ops, func_name='suite2p_registration', file_name='ops')

    return info
