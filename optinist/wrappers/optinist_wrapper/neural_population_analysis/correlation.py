from optinist.wrappers.data_wrapper import *
from optinist.api.nwb.const import NWBDATASET


def correlation(
        neural_data: FluoData,
        iscell: IscellData=None,
        nwbfile: NWBFile=None,
        params: dict=None
    ) -> dict():

    neural_data = neural_data.data

    # data shold be time x component matrix
    if params['transpose']:
        X = neural_data.transpose()
    else:
        X = neural_data

    if iscell is not None:
        iscell = iscell.data
        ind  = np.where(iscell > 0)[0]
        X = X[ind, :]

    num_cell = X.shape[0]

    # calculate correlation
    corr = np.corrcoef(X)
    for i in range(num_cell):
        corr[i, i] = np.nan

    info = {}
    info['corr'] = CorrelationData(
        corr,
        file_name='corr'
    )

    # NWB追加
    if nwbfile is not None:
        nwbfile[NWBDATASET.POSTPROCESS] = {
            'corr': corr,
        }

    info['nwbfile'] = nwbfile

    return info
