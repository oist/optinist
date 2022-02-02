from wrappers.data_wrapper import *
from wrappers.args_check import args_check
from wrappers.optinist_wrapper.utils import standard_norm

def PCA(
        neural_data: TimeSeriesData,
        iscell: IscellData=None,
        nwbfile: NWBFile=None,
        params: dict=None
    ) -> {}:
    # modules specific to function
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA

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

    # # preprocessing  ##################
    tX = standard_norm(X, params['standard_mean'], params['standard_std'])

    # calculate PCA ##################
    pca = PCA(**params['PCA'])
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
