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
        timeseries1: TimeSeriesData,
        behaviors: TimeSeriesData, 
        iscell: IscellData=None,
        params: dict=None
    ) -> {'actual_predicted': TimeSeriesData}:

    # modules specific to function
    import statsmodels.api as sm
    from sklearn.preprocessing import StandardScaler
    import pandas as pd

    timeseries1 = timeseries1.data
    behaviors = behaviors.data

    if iscell is not None:
        iscell = iscell.data
        ind  = np.where(iscell > 0)[0]
        timeseries1 = timeseries1[ind, :]
        behaviors = behaviors[ind, :]

    # data shold be time x component matrix
    neural_data = timeseries1 # neural data
    behavioral_data = behaviors[:, params['target_index']].reshape(-1, 1)

    # preprocessing
    neural_data = standard_norm(neural_data, params['standard_n_mean'], params['standard_n_std'])
    behavioral_data = standard_norm(behavioral_data, params['standard_b_mean'], params['standard_b_std'])

    # calculate GLM
    neural_data = pd.DataFrame(neural_data)
    behavioral_data = pd.DataFrame(behavioral_data)

    if(params['add_constant']):
        neural_data = sm.add_constant(neural_data, prepend=False)

    # set link function
    link = getattr(sm.genmod.families.links, params['link'])()

    # set family
    family = eval('sm.families.' + params['family'] + '(link=link)')

    # model fit
    model = sm.GLM(
        behavioral_data, neural_data, family=family, offset=params['offset'], 
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
