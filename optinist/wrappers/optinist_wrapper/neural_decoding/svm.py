from wrappers.data_wrapper import *
from wrappers.args_check import args_check
from wrappers.optinist_wrapper.utils import standard_norm

def SVM(
        neural_data: TimeSeriesData,
        behaviors_data: TimeSeriesData,
        iscell: IscellData=None,
        nwbfile: NWBFile=None,
        params: dict=None
    ) -> {}:

    # modules specific to function
    import numpy as np
    from sklearn import svm
    from sklearn.model_selection import GridSearchCV
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import StratifiedKFold

    neural_data = neural_data.data
    behaviors_data = behaviors_data.data

    # data shold be time x component matrix
    if params['transpose_x']:
        X = neural_data.transpose()
    else:
        X = neural_data

    if params['transpose_y']:
        Y = behaviors_data.transpose()
    else:
        Y = behaviors_data

    assert X.shape[0] == Y.shape[0], f"""
        neural_data and behaviors_data is not same dimension,
        neural.shape{X.shape}, behavior.shape{Y.shape}"""

    if iscell is not None:
        iscell = iscell.data
        ind  = np.where(iscell > 0)[0]
        X = X[ind, :]
        Y = Y[ind, :]

    Y = Y[:, params['target_index']].reshape(-1, 1)

    # # preprocessing  ##################
    tX = standard_norm(X, params['standard_x_mean'], params['standard_x_std'])

    hp = params['SVC']

    # SVM determination of hyper parameters if needed ##################
    gs_clf=[]
    if params['use_grid_search']:
        param_grid = [params['grid_search']['param_grid']]
        gs_clf = GridSearchCV(svm.SVC(), param_grid=param_grid, **params['grid_search']['CV'])

        gs_clf.fit(tX, Y)

        # insert best parameter to hp dictionary
        keys = list(gs_clf.best_params_.keys())
        for i in range(len(keys)):
            hp[keys[i]] = gs_clf.best_params_[keys[i]]
        gs_best_params = gs_clf.best_params_

    # cross validation of SVM using best grid search paraneters ##################
    skf = StratifiedKFold(**params['CV'])

    score = []
    classifier =[]
    for train_index, test_index in skf.split(tX, Y):

        clf = svm.SVC(**hp)

        if tX.shape[0] == 1:
            clf.fit(tX[train_index].reshape(-1, 1), Y[train_index])
            score.append(clf.score(tX[test_index].reshape(-1, 1), Y[test_index]))
        else:
            clf.fit(tX[train_index, :], Y[train_index])
            score.append(clf.score(tX[test_index, :], Y[test_index]))

        classifier.append(clf)

    info = {}
    info['score'] = BarData(
        score,
        func_name='svm',
        file_name='score'
    )

    # NWB追加
    if nwbfile is not None:
        nwbfile['add_postprocess'] = {
            'score': score,
        }

    info['nwbfile'] = nwbfile

    return info
