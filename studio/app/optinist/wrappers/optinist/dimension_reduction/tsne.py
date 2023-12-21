from studio.app.common.dataclass import ScatterData
from studio.app.optinist.core.nwb.nwb import NWBDATASET
from studio.app.optinist.dataclass import FluoData, IscellData
from studio.app.optinist.wrappers.optinist.utils import standard_norm


def TSNE(
    neural_data: FluoData,
    output_dir: str,
    iscell: IscellData = None,
    params: dict = None,
    **kwargs,
) -> dict():
    import numpy as np
    from sklearn.manifold import TSNE

    function_id = output_dir.split("/")[-1]
    print("start TSNE:", function_id)

    neural_data = neural_data.data

    # data should be time x component matrix
    if params["transpose"]:
        X = neural_data.transpose()
    else:
        X = neural_data

    if iscell is not None:
        iscell = iscell.data
        ind = np.where(iscell > 0)[0]
        X = X[:, ind]

    # preprocessing
    tX = standard_norm(X, params["standard_mean"], params["standard_std"])

    # calculate TSNE
    tsne = TSNE(**params["TSNE"])

    proj_X = tsne.fit_transform(tX)

    # NWB追加
    nwbfile = {}
    nwbfile[NWBDATASET.POSTPROCESS] = {function_id: {"projectedNd": proj_X}}

    info = {
        "projectedNd": ScatterData(proj_X, file_name="projectedNd"),
        "nwbfile": nwbfile,
    }

    return info
