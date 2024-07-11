from studio.app.common.core.experiment.experiment import ExptOutputPathIds
from studio.app.common.core.logger import AppLogger
from studio.app.common.dataclass import TimeSeriesData
from studio.app.optinist.core.nwb.nwb import NWBDATASET
from studio.app.optinist.dataclass import FluoData, IscellData

logger = AppLogger.get_logger()


def cross_correlation(
    neural_data: FluoData,
    output_dir: str,
    iscell: IscellData = None,
    params: dict = None,
    **kwargs,
) -> dict():
    import itertools

    import numpy as np
    import scipy.signal as ss
    import scipy.stats as stats
    from tqdm import tqdm

    function_id = ExptOutputPathIds(output_dir).function_id
    logger.info("start cross_correlation: %s", function_id)

    neural_data = neural_data.data

    # data should be time x component matrix
    if params["transpose"]:
        X = neural_data.transpose()
    else:
        X = neural_data

    if iscell is not None:
        iscell = iscell.data
        ind = np.where(iscell > 0)[0]
        X = X[ind, :]

    # calculate cross correlation
    num_cell = X.shape[0]
    data_len = X.shape[1]
    shuffle_num = params["shuffle_sample_number"]

    lags = ss.correlation_lags(data_len, data_len, mode="same")
    ind = np.intersect1d(
        np.where(-params["lags"] <= lags), np.where(lags <= params["lags"])
    )
    x = lags[ind]

    mat = np.zeros([num_cell, num_cell, len(x)])
    s_confint = np.zeros([num_cell, num_cell, len(x), 2])
    s_mean = np.zeros([num_cell, num_cell, len(x)])

    for i in tqdm(range(num_cell)):
        # for i in tqdm(range(3)):
        for j in tqdm(range(num_cell)):
            ccvals = ss.correlate(
                X[i, :], X[j, :], method=params["method"], mode="same"
            )
            mat[i, j, :] = ccvals[ind]

            # baseline
            tp_s_mat = np.zeros([len(x), shuffle_num])
            for k in range(shuffle_num):
                tp = X[j, :]
                np.random.shuffle(tp)
                ccvals = ss.correlate(X[i, :], tp, method=params["method"], mode="same")
                tp_s_mat[:, k] = ccvals[ind]

            b_mean = np.mean(tp_s_mat, axis=1)
            b_sem = stats.sem(tp_s_mat, axis=1)
            b_df = shuffle_num - 1

            interval = np.zeros([len(x), 2])
            for k in range(len(x)):
                interval[k, :] = stats.t.interval(
                    params["shuffle_confidence_interval"], b_df, loc=0, scale=b_sem[k]
                )

            s_confint[i, j, :, 0] = interval[:, 0]
            s_confint[i, j, :, 1] = interval[:, 1]
            s_mean[i, j, :] = b_mean

    # NWB追加
    nwbfile = {}
    nwbfile[NWBDATASET.POSTPROCESS] = {
        function_id: {
            "mat": mat,
            "baseline": s_mean,
            "base_confint": s_confint,
        }
    }

    info = {
        # 'cross_correlation': ScatterData(mat),
        "nwbfile": nwbfile
    }

    # output structures
    cb = list(itertools.combinations(range(num_cell), 2))

    for i in range(len(cb)):
        arr1 = np.stack([x, mat[cb[i][0], cb[i][1], :]], axis=1)
        arr2 = np.stack(
            [
                x,
                s_mean[cb[i][0], cb[i][1], :],
                s_confint[cb[i][0], cb[i][1], :, 0],
                s_confint[cb[i][0], cb[i][1], :, 1],
            ],
            axis=1,
        )

    name = f"{str(cb[i][0])}-{str(cb[i][1])}"
    info[name] = TimeSeriesData(arr1.T, file_name=name)
    name = f"shuffle {str(cb[i][0])}-{str(cb[i][1])}"
    info[name] = TimeSeriesData(arr2.T, file_name=name)

    return info
