from wrappers.data_wrapper import *
from wrappers.args_check import args_check
from wrappers.optinist_wrapper.utils import standard_norm


@args_check
def PCA(
        timeseries: TimeSeriesData, iscell: IscellData=None, params: dict=None
    ) -> {}:
    # modules specific to function
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA

    timeseries = timeseries.data
    

    if iscell is not None:
        iscell = iscell.data
        ind = np.where(iscell > 0)[0]
        timeseries = timeseries[ind, :]

    # data shold be time x component matrix
    X = timeseries.transpose()

    # # preprocessing  ##################
    tX = standard_norm(X, params['standard_mean'], params['standard_std'])

    # calculate PCA ##################
    pca = PCA(
        n_components=params['n_components'],
        copy=params['copy'],
        whiten= params['whiten'],
        svd_solver= params['svd_solver'],
        tol= params['tol'],
        iterated_power=params['iterated_power'],
        random_state= params['random_state']
    )
    proj_X = pca.fit_transform(tX)

    info = {}
    info['components'] = CorrelationData(
        pca.components_,
        func_name='pca',
        file_name='components'
    )
    info['explained_variance_ratio'] = TimeSeriesData(
        pca.explained_variance_ratio_,
        func_name='pca',
        file_name='evr'
    )
    info['projected2d'] = ScatterData(
        proj_X[:, 0:2],
        func_name='pca',
        file_name='projected2d'
    )  # change to 2D scatter plots

    return info
