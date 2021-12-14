from wrappers.data_wrapper import *
from wrappers.args_check import args_check


@args_check
def suite2p_spike_deconv(ops: Suite2pData, params: dict=None) -> {'ops': Suite2pData, 'spks': TimeSeriesData}:
    ops = ops.data
    from suite2p import extraction

    ops = {**ops, **params}

    dF = ops['F'].copy() - ops['neucoeff']*ops['Fneu']
    dF = extraction.preprocess(
        F=dF,
        baseline=ops['baseline'],
        win_baseline=ops['win_baseline'],
        sig_baseline=ops['sig_baseline'],
        fs=ops['fs'],
        prctile_baseline=ops['prctile_baseline']
    )
    spks = extraction.oasis(F=dF, batch_size=ops['batch_size'], tau=ops['tau'], fs=ops['fs'])

    ops['spks'] = spks

    info = {}
    info['ops'] = Suite2pData(ops)
    info['spks'] = TimeSeriesData(spks, func_name='sute2p_spike_deconv', file_name='spks')

    return info
