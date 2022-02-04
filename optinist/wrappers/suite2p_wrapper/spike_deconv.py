from wrappers.data_wrapper import *
from wrappers.args_check import args_check


def suite2p_spike_deconv(
        ops: Suite2pData, nwbfile: NWBFile=None, params: dict=None
    ) -> {'ops': Suite2pData, 'spks': TimeSeriesData}:
    from suite2p import extraction, default_ops
    print('start suite2_spike_deconv')
    ops = ops.data

    ops = {**default_ops(), **ops, **params}

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

    # NWBを追加
    ### Fluorenceを追加
    if nwbfile is not None:
        nwbfile['add_fluorescence'] = {}
        for name, data in zip(['Fluorescence', 'Neuropil', 'Deconvolved'], [ops['F'], ops['Fneu'], spks]):
            nwbfile['add_fluorescence'][name] = {
                'table_name': name,
                'region': list(range(len(data))),
                'name': name,
                'data': data,
                'unit': 'lumens',
                'rate': ops['fs'],
            }

    info = {}
    info['ops'] = Suite2pData(ops)
    info['spks'] = TimeSeriesData(spks, func_name='sute2p_spike_deconv', file_name='spks')
    info['nwbfile'] = nwbfile

    return info
