from wrappers.data_wrapper import *
from wrappers.args_check import args_check
from wrappers.optinist_wrapper.utils import standard_norm


@args_check
def CCA(
        timeseries: TimeSeriesData,
        behaviors: TimeSeriesData,
        iscell: IscellData=None,
        params: dict=None
    ) -> {}:
    from sklearn.cross_decomposition import CCA

    timeseries = timeseries.data
    behaviors = behaviors.data

    if iscell is not None:
        iscell = iscell.data
        ind = np.where(iscell > 0)[0]
        timeseries = timeseries[ind, :]
        behaviors = behaviors[ind, :]

    # data shold be time x component matrix
    if params['transpose_x']:
        X = timeseries.transpose()
    else:
        X = timeseries

    if params['transpose_y']:
        Y = behaviors.transpose()
    else:
        Y = behaviors

    Y = Y[:, params['target_index']].reshape(-1, 1)

    assert X.shape[0] == Y.shape[0], f'X and Y is not same data, X.shape{X.shape}, Y.shape{Y.shape}'

    # # preprocessing  ##################
    tX = standard_norm(X, params['standard_x_mean'], params['standard_x_std'])
    tY = standard_norm(Y, params['standard_y_mean'], params['standard_y_std'])

    # calculate CCA 
    cca = CCA(**params['CCA'])
    projX, projY = cca.fit_transform(tX, tY)

    proj = np.concatenate([projX, projY], axis=1)

    info = {}
    info['projected2d'] = ScatterData(
        proj,
        func_name='cca',
        file_name='projected2d'
    )

    return info
