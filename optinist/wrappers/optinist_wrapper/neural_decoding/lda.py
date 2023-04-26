from optinist.api.dataclass.dataclass import BarData, BehaviorData, FluoData, IscellData
from optinist.api.nwb.nwb import NWBDATASET
from optinist.wrappers.optinist_wrapper.utils import standard_norm


def LDA(
    neural_data: FluoData,
    behaviors_data: BehaviorData,
    output_dir: str,
    iscell: IscellData = None,
    params: dict = None,
) -> dict():
    # modules specific to function
    import numpy as np
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
    from sklearn.model_selection import StratifiedKFold

    neural_data = neural_data.data
    behaviors_data = behaviors_data.data

    # data shold be time x component matrix
    if params["transpose_x"]:
        X = neural_data.transpose()
    else:
        X = neural_data

    if params["transpose_y"]:
        Y = behaviors_data.transpose()
    else:
        Y = behaviors_data

    assert (
        X.shape[0] == Y.shape[0]
    ), f"""
        neural_data and behaviors_data is not same dimension,
        neural.shape{neural_data.shape}, behavior.shape{behaviors_data.shape}"""

    if iscell is not None:
        iscell = iscell.data
        ind = np.where(iscell > 0)[0]
        X = X[:, ind]

    Y = Y[:, params["target_index"]].reshape(-1, 1)

    # preprocessing
    tX = standard_norm(X, params["standard_x_mean"], params["standard_x_std"])

    # cross validation of LDA model
    skf = StratifiedKFold(**params["CV"])

    score = []
    classifier = []
    for train_index, test_index in skf.split(tX, Y):
        clf = LDA(**params["LDA"])

        if tX.shape[0] == 1:
            clf.fit(tX[train_index].reshape(-1, 1), Y[train_index])
            score.append(clf.score(tX[test_index].reshape(-1, 1), Y[test_index]))
            classifier.append(clf)
        else:
            clf.fit(tX[train_index, :], Y[train_index])
            score.append(clf.score(tX[test_index, :], Y[test_index]))
            classifier.append(clf)

    # NWB追加
    nwbfile = {}
    nwbfile[NWBDATASET.POSTPROCESS] = {
        "score": score,
    }

    info = {"score": BarData(score, file_name="score"), "nwbfile": nwbfile}

    return info
