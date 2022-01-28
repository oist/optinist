from wrappers.data_wrapper import *
from wrappers.args_check import args_check
from wrappers.optinist_wrapper.utils import standard_norm

@args_check
def TSNE(
        timeseries: TimeSeriesData,
        iscell: IscellData=None,
        params: dict=None
    ) -> {}:

    from sklearn.manifold import TSNE

    timeseries = timeseries.data

    if iscell is not None:
        iscell = iscell.data
        ind = np.where(iscell > 0)[0]
        timeseries = timeseries[ind, :]

    # data shold be time x component matrix
    if params['transpose']:
        X = timeseries.transpose()
    else:
        X = timeseries

    # preprocessing  ##################
    tX = standard_norm(X, params['standard_mean'], params['standard_std'])

    # calculate TSNE ##################
    tsne = TSNE(**params['TSNE'])

    proj_X = tsne.fit_transform(tX)

    info = {}
    info['projected2d'] = ScatterData(
        proj_X[:, 0:2],
        func_name='tsne',
        file_name='projected2d'
    )

    return info
