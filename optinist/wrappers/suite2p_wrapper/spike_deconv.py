from optinist.api.dataclass.dataclass import *
from optinist.api.nwb.nwb import NWBDATASET


def suite2p_spike_deconv(
        ops: Suite2pData,
        nwbfile: NWBFile=None,
        params: dict=None
    ) -> dict(ops=Suite2pData, spks=FluoData):
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
        if NWBDATASET.FLUORESCENCE not in nwbfile.keys():
            nwbfile[NWBDATASET.FLUORESCENCE] = {}
        name = 'Deconvolved'
        nwbfile[NWBDATASET.FLUORESCENCE][name] = {
            'table_name': name,
            'region': list(range(len(spks))),
            'name': name,
            'data': spks,
            'unit': 'lumens',
            'rate': ops['fs'],
        }

    info = {
        'ops': Suite2pData(ops),
        'spks': FluoData(spks, file_name='spks'),
        'nwbfile': nwbfile
    }

    return info
