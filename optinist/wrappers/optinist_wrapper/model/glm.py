#  decoding neural activity by GLM
#  input:  A:matrix[num_cell x timeseries ]   B:timeseries(behavior)[1 x timeseries]
#  Generalized linear modeling using statsmodels
#
#　https://www.statsmodels.org/stable/glm.html

from wrappers.data_wrapper import *
from wrappers.args_check import args_check
from wrappers.optinist_wrapper.utils import standard_norm

@args_check
def GLM(
        neural_data: TimeSeriesData,
        behaviors_data: TimeSeriesData, 
        iscell: IscellData=None,
        params: dict=None
    ) -> {'actual_predicted': TimeSeriesData}:

    # modules specific to function
    import statsmodels.api as sm
    from sklearn.preprocessing import StandardScaler
    import pandas as pd

    neural_data = neural_data.data
    behaviors_data = behaviors_data.data

    if iscell is not None:
        iscell = iscell.data
        ind  = np.where(iscell > 0)[0]
        neural_data = neural_data[ind, :]
        behaviors_data = behaviors_data[ind, :]

    # data shold be time x component matrix
    if params['transpose_x']:
        X = neural_data.transpose()
    else:
        X = neural_data

    if params['transpose_y']:
        Y = behaviors_data.transpose()
    else:
        Y = behaviors_data

    Y = Y[:, params['target_index']].reshape(-1, 1)

    assert X.shape[0] == Y.shape[0], f'X and Y is not same data, X.shape{X.shape}, Y.shape{Y.shape}'

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
    model = sm.GLM(
        tY, tX, family=family, offset=params['offset'], 
        exposure=params['exposure'], missing=params['missing'])
    Res = model.fit()

    # main results for plot
    # plot should be reconsidered --- what they should be!
    info = {}
    #info['params'] = BarSeriesData(Res.params.values)  # add something 1d but not timesereies
    info['actual_predicted'] = ScatterData(
        np.array([Res._endog, Res.mu]).transpose(),
        func_name='glm',
        file_name='actual_predicted'
    )

    return info
