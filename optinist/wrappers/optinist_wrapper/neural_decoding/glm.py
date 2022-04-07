#  decoding neural activity by GLM
#  input:  A:matrix[num_cell x timeseries ]   B:timeseries(behavior)[1 x timeseries]
#  Generalized linear modeling using statsmodels
#
#　https://www.statsmodels.org/stable/glm.html

from optinist.wrappers.data_wrapper import *
from optinist.wrappers.optinist_wrapper.utils import standard_norm
from optinist.wrappers.nwb_wrapper.const import NWBDATASET

def GLM(
        neural_data: FluoData,
        behaviors_data: BehaviorData,
        iscell: IscellData=None,
        nwbfile: NWBFile=None,
        params: dict=None
    ) -> {}:

    # modules specific to function
    import statsmodels.api as sm
    from sklearn.preprocessing import StandardScaler
    import pandas as pd

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
        X = X[:, ind]
        Y = Y[:, ind]

    Y = Y[:, params['target_index']].reshape(-1, 1)

    # preprocessing
    tX = standard_norm(X, params['standard_x_mean'], params['standard_x_std'])
    tY = standard_norm(Y, params['standard_y_mean'], params['standard_y_std'])

    # calculate GLM
    tX = pd.DataFrame(tX)
    tY = pd.DataFrame(tY)

    if(params['add_constant']):
        tX = sm.add_constant(tX, prepend=False)

    # set link function
    link = getattr(sm.genmod.families.links, params['link'])()

    # set family
    family = eval('sm.families.' + params['family'] + '(link=link)')

    # model fit
    model = sm.GLM(tY, tX, family=family, **params['GLM'])
    Res = model.fit()

    # main results for plot
    # plot should be reconsidered --- what they should be!
    info = {}
    #info['params'] = BarSeriesData(Res.params.values)  # add something 1d but not timesereies
    info['actual_predicted'] = ScatterData(
        np.array([Res._endog, Res.mu]).transpose(),
        file_name='actual_predicted'
    )
    info['params'] = BarData(
        Res.params.values,
        file_name='params'
    )
    info['summary'] = HTMLData(
        Res.summary().as_html(),
        file_name='summary',
    )

    # NWB追加
    if nwbfile is not None:
        nwbfile[NWBDATASET.POSTPROCESS] = {
            'actual_predicted': np.array([Res._endog, Res.mu]).transpose(),
            'params': Res.params.values,
        }

    info['nwbfile'] = nwbfile

    return info
