from wrappers.data_wrapper import *
from wrappers.args_check import args_check

@args_check
def correlation(
        timeseries: TimeSeriesData, iscell: IscellData=None, params: dict=None
    ) -> {'corr': CorrelationData}:

    timeseries = timeseries.data

    if iscell is not None:
        iscell = iscell.data
        ind  = np.where(iscell > 0)[0]
        timeseries = timeseries[ind, :]

    # data shold be time x component matrix
    if params['transpose']:
        X = timeseries.transpose()
    else:
        X = timeseries

    num_cell = X.shape[0]

    # calculate correlation
    corr = np.corrcoef(X)
    for i in range(num_cell):
        corr[i, i] = np.nan

    info = {}
    info['corr'] = CorrelationData(
        corr,
        func_name='correlation',
        file_name='corr'
    )

    return info
