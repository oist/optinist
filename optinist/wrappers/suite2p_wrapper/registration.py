def suite2p_registration(ops, opts=None):
    import numpy as np
    from suite2p import registration

    ######### REGISTRATION #########
    refImg = ops['refImg'] if 'refImg' in ops and ops.get('force_refImg', False) else None
    ops = registration.register_binary(ops, refImg=refImg) # register binary

    # compute metrics for registration
    if ops.get('do_regmetrics', True) and ops['nframes']>=1500:
        ops = registration.get_pc_metrics(ops)

    ops['images'] = ops['refImg'].astype(np.uint8)
    ops['meanImgE'] = ops['meanImgE']

    return ops
