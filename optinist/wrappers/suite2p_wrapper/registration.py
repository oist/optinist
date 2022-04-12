from optinist.wrappers.data_wrapper import *


def suite2p_registration(
        ops: Suite2pData,
        nwbfile: NWBFile=None,
        params: dict=None
    ) -> dict(ops=Suite2pData):
    import numpy as np
    from suite2p import registration, io, default_ops, io
    ops = ops.data
    refImg = ops['meanImg']
    print('start suite2_registration')

    ######### REGISTRATION #########
    if len(refImg.shape) == 3:
        refImg = refImg[0]

    ops = {**default_ops(), **ops, **params}

    # register binary
    ops = registration.register_binary(ops, refImg=refImg)

    # compute metrics for registration
    if ops.get('do_regmetrics', True) and ops['nframes']>=1500:
        ops = registration.get_pc_metrics(ops)

    info = {}
    info['refImg'] = ImageData(ops['refImg'], file_name='refImg')
    info['meanImgE'] = ImageData(ops['meanImgE'], file_name='meanImgE')
    info['ops'] = Suite2pData(ops, file_name='ops')

    return info
