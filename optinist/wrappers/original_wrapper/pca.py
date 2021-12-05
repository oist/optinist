from wrappers.data_wrapper import *
from wrappers.args_check import args_check

@args_check
def PCA(timeseries: TimeSeriesData, iscell: IscellData, params: dict=None):
    # modules specific to function
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA

    timeseries = timeseries.data
    iscell = iscell.data
    ind = np.where(iscell > 0)[0]
    timeseries = timeseries[ind, :]

    # data shold be time x component matrix
    X = timeseries.transpose()

    # preprocessing  ##################
    sc = StandardScaler(with_mean = params['standardization_mean'], with_std = params['standardization_std'])
    tX = sc.fit_transform(X)

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
    pca.fit(tX)

    proj_X = pca.transform(tX)

    # output structures ###################
    Out_nwb = {}
    Out_nwb['components'] = pca.components_
    Out_nwb['explained_variance'] = pca.explained_variance_
    Out_nwb['explained_variance_ratio'] = pca.explained_variance_ratio_
    Out_nwb['singular_values'] = pca.singular_values_
    Out_nwb['mean'] = pca.mean_
    Out_nwb['n_components'] = pca.n_components_
    Out_nwb['n_samples'] = pca.n_samples_
    Out_nwb['noise_variance'] = pca.noise_variance_
    Out_nwb['n_features_in'] = pca.n_features_in_
    #Out_nwb['feature_names_in'] = pca.feature_names_in_

    info = {}
    info['components'] = CorrelationData(pca.components_, 'components')
    info['explained_variance_ratio'] = TimeSeriesData(pca.explained_variance_ratio_, 'evr')
    info['projected'] = TimeSeriesData(proj_X[:, 0:2].transpose(1, 0), 'projected')  # change to 2D scatter plots

    return info
