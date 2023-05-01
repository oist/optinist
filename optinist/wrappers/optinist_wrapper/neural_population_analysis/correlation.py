from optinist.api.dataclass.dataclass import FluoData, HeatMapData, IscellData
from optinist.api.nwb.nwb import NWBDATASET


def correlation(
    neural_data: FluoData,
    output_dir: str,
    iscell: IscellData = None,
    params: dict = None,
) -> dict():
    import numpy as np

    neural_data = neural_data.data

    # data shold be time x component matrix
    if params["transpose"]:
        X = neural_data.transpose()
    else:
        X = neural_data

    if iscell is not None:
        iscell = iscell.data
        ind = np.where(iscell > 0)[0]
        X = X[ind, :]

    num_cell = X.shape[0]

    # calculate correlation
    corr = np.corrcoef(X)
    for i in range(num_cell):
        corr[i, i] = np.nan

    # NWB追加
    nwbfile = {}
    nwbfile[NWBDATASET.POSTPROCESS] = {
        "corr": corr,
    }

    info = {
        "corr": HeatMapData(corr, file_name="corr"),
        "nwbfile": nwbfile,
    }

    return info
