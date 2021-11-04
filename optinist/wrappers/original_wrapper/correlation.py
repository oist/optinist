from wrappers.data_wrapper import *
from wrappers.args_check import args_check

@args_check
def correlation(timeseries: TimeSeriesData, iscell: IscellData):
    timeseries = timeseries.data
    iscell = iscell.data

    ind  = np.where(iscell > 0)[0]

    timeseries = timeseries[ind, :]

    # calculate correlation ##################
    num_cell = timeseries.shape[0]

    corr = np.ones([num_cell, num_cell]) * np.nan
    for i in range(num_cell):
        for j in range(num_cell):
            if(i != j):
                corr[i, j] = np.corrcoef(timeseries[i, :], timeseries[j, :])[0, 1]

    info = {}
    info['corr'] = CorrelationData(corr)

    return info
