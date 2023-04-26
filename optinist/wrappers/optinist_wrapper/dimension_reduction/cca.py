from optinist.api.dataclass.dataclass import (
    BarData,
    BehaviorData,
    FluoData,
    IscellData,
    ScatterData,
)
from optinist.api.nwb.nwb import NWBDATASET
from optinist.wrappers.optinist_wrapper.utils import standard_norm


def CCA(
    neural_data: FluoData,
    behaviors_data: BehaviorData,
    output_dir: str,
    iscell: IscellData = None,
    params: dict = None,
) -> dict():
    import numpy as np
    from sklearn.cross_decomposition import CCA

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
        neural.shape{X.shape}, behavior.shape{Y.shape}"""

    if iscell is not None:
        iscell = iscell.data
        ind = np.where(iscell > 0)[0]
        X = X[:, ind]

    Y = Y[:, params["target_index"]].reshape(-1, 1)

    # preprocessing
    tX = standard_norm(X, params["standard_x_mean"], params["standard_x_std"])
    tY = standard_norm(Y, params["standard_y_mean"], params["standard_y_std"])

    # calculate CCA
    cca = CCA(**params["CCA"])
    projX, projY = cca.fit_transform(tX, tY)

    proj = np.concatenate([projX, projY], axis=1)

    # NWB追加
    nwbfile = {}
    nwbfile[NWBDATASET.POSTPROCESS] = {
        "projectedNd": proj,
        "x_weights": cca.x_weights_,  # singular vectors
        "y_weights": cca.y_weights_,
        "x_loadings_": cca.x_rotations_,
        "y_loadings_": cca.x_rotations_,
        "coef": cca.coef_,
        "n_iter_": cca.n_iter_,
        # 'n_features_in_': [cca.n_features_in_],
    }

    info = {
        "projectedNd": ScatterData(proj, file_name="projectedNd"),
        "coef": BarData(cca.coef_.flatten(), file_name="coef"),
        "nwbfile": nwbfile,
    }

    return info
