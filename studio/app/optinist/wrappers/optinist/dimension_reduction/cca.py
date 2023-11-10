from studio.app.common.dataclass import BarData, ScatterData
from studio.app.optinist.core.nwb.nwb import NWBDATASET
from studio.app.optinist.dataclass import BehaviorData, FluoData, IscellData
from studio.app.optinist.wrappers.optinist.utils import standard_norm


def CCA(
    neural_data: FluoData,
    behaviors_data: BehaviorData,
    output_dir: str,
    iscell: IscellData = None,
    params: dict = None,
    **kwargs,
) -> dict():
    import numpy as np
    from sklearn.cross_decomposition import CCA

    function_id = output_dir.split("/")[-1]
    print("start cca:", function_id)

    neural_data = neural_data.data
    behaviors_data = behaviors_data.data

    # data should be time x component matrix
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
        function_id: {
            "projectedNd": proj,
            "x_weights": cca.x_weights_,  # singular vectors
            "y_weights": cca.y_weights_,
            "x_loadings_": cca.x_rotations_,
            "y_loadings_": cca.x_rotations_,
            "coef": cca.coef_,
            "n_iter_": cca.n_iter_,
            # 'n_features_in_': [cca.n_features_in_],
        }
    }

    info = {
        "projectedNd": ScatterData(proj, file_name="projectedNd"),
        "coef": BarData(cca.coef_.flatten(), file_name="coef"),
        "nwbfile": nwbfile,
    }

    return info
