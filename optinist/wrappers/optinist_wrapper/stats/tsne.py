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

    # calculate PCA ##################
    tsne = TSNE(
        n_components=params['n_components'],
        perplexity=params['perplexity'],
        early_exaggeration= params['early_exaggeration'],
        learning_rate= params['learning_rate'],
        n_iter=params['n_iter'],
        n_iter_without_progress= params['n_iter_without_progress'],
        min_grad_norm= params['min_grad_norm'],
        metric=params['metric'],
        init=params['init'],
        random_state=params['random_state'],
        method=params['method'],
        angle=params['angle'],
        n_jobs=params['n_jobs'],
        square_distances=params['square_distances']
    )

    proj_X = tsne.fit_transform(tX)

    info = {}
    info['projected2d'] = ScatterData(
        proj_X[:, 0:2],
        func_name='tsne',
        file_name='projected2d'
    )

    return info
