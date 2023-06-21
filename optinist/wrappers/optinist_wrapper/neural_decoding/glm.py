#  decoding neural activity by GLM
#  input:  A:matrix[num_cell x timeseries ]   B:timeseries(behavior)[1 x timeseries]
#  Generalized linear modeling using statsmodels
#
#  https://www.statsmodels.org/stable/glm.html

from optinist.api.dataclass.dataclass import (
    BarData,
    BehaviorData,
    FluoData,
    HTMLData,
    IscellData,
    ScatterData,
)
from optinist.api.nwb.nwb import NWBDATASET
from optinist.wrappers.optinist_wrapper.utils import standard_norm


def GLM(
    neural_data: FluoData,
    behaviors_data: BehaviorData,
    output_dir: str,
    iscell: IscellData = None,
    params: dict = None,
) -> dict():
    # modules specific to function
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm

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
        Y = Y[:, ind]

    Y = Y[:, params["target_index"]].reshape(-1, 1)

    # preprocessing
    tX = standard_norm(X, params["standard_x_mean"], params["standard_x_std"])
    tY = standard_norm(Y, params["standard_y_mean"], params["standard_y_std"])

    # calculate GLM
    tX = pd.DataFrame(tX)
    tY = pd.DataFrame(tY)

    if params["add_constant"]:
        tX = sm.add_constant(tX, prepend=False)

    # set family
    link = getattr(sm.genmod.families.links, params["link"])()  # noqa: F841
    family = eval(f"sm.families.{params['family']}(link=link)")

    # model fit
    model = sm.GLM(tY, tX, family=family, **params["GLM"])
    Res = model.fit()

    # NWB追加
    nwbfile = {}
    nwbfile[NWBDATASET.POSTPROCESS] = {
        "actual_predicted": np.array([Res._endog, Res.mu]).transpose(),
        "params": Res.params.values,
        "pvalues": Res.pvalues.values,
        "tvalues": Res.tvalues.values,  # z
        "aic": [Res.aic],
        "bic_llf": [Res.bic_llf],
        "llf": [Res.llf],  # log-Likelihood
        "pearson_chi2": [Res.pearson_chi2],
        "df_model": [Res.df_model],
        "df_resid": [Res.df_resid],
        "scale": [Res.scale],
        "mu": Res.mu,
        "endog": Res._endog,
    }

    # main results for plot
    # plot should be reconsidered --- what they should be!
    info = {
        "actual_predicted": ScatterData(
            np.array([Res._endog, Res.mu]).transpose(), file_name="actual_predicted"
        ),
        "params": BarData(Res.params.values, file_name="params"),
        "textout": HTMLData(Res.summary().as_html(), file_name="textout"),
        "nwbfile": nwbfile,
    }

    return info
