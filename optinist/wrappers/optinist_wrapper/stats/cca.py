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
    X = timeseries.transpose()
    Y = behaviors.transpose()

    # # preprocessing  ##################
    tX = standard_norm(X, params['standard_x_mean'], params['standard_x_std'])
    tY = standard_norm(Y, params['standard_y_mean'], params['standard_y_std'])

    # calculate CCA 
    cca = CCA(
        n_components=params['n_components'],
        scale=params['scale'],
        max_iter=params['max_iter'],
        tol=params['tol'],
        copy=params['copy']
    )
    projX, projY = cca.fit_transform(tX, tY)

    info = {}
    info['projected2d'] = ScatterData(
        proj_X[:, 0:2],
        func_name='cca',
        file_name='projected2d'
    )

    return info
