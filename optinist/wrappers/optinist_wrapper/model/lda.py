from wrappers.data_wrapper import *
from wrappers.args_check import args_check
from wrappers.optinist_wrapper.utils import standard_norm


@args_check
def LDA(
        timeseries: TimeSeriesData,
        behaviors: TimeSeriesData,
        iscell: IscellData=None,
        params: dict=None
    ):

    # modules specific to function
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import StratifiedKFold

    timeseries = timeseries.data
    behaviors = behaviors.data

    if iscell is not None:
        iscell = iscell.data
        ind = np.where(iscell > 0)[0]
        timeseries = timeseries[ind, :]
        behaviors = behaviors[ind, :]

    # # preprocessing  ##################
    if params['transpose_x']:
        X = timeseries.transpose()
    else:
        X = timeseries

    if params['transpose_y']:
        Y = behaviors.transpose()
    else:
        Y = behaviors

    Y = Y[:, params['target_index']].reshape(-1, 1)

    # # preprocessing  ##################
    tX = standard_norm(X, params['standard_x_mean'], params['standard_x_std'])

    # cross validation of LDA model ##################
    skf = StratifiedKFold(**params['CV'])

    score = []
    classifier = []
    for train_index, test_index in skf.split(tX, Y):
        clf = LDA(**params['LDA'])


        if (tX.shape[0] == 1):
            clf.fit(tX[train_index].reshape(-1, 1), Y[train_index])
            score.append(clf.score(tX[test_index].reshape(-1, 1), Y[test_index]))
            classifier.append(clf)
        else:
            clf.fit(tX[train_index, :], Y[train_index])
            score.append(clf.score(tX[test_index, :], Y[test_index]))
            classifier.append(clf)

    info = {}
    info['score'] = score

    return info
