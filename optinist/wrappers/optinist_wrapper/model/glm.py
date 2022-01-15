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
        timeseries2: TimeSeriesData, 
        iscell: IscellData=None,
        params: dict=None
    ) -> {'actual_predicted': TimeSeriesData}:

    # modules specific to function
    import statsmodels.api as sm
    from sklearn.preprocessing import StandardScaler
    import pandas as pd

    timeseries1 = timeseries1.data
    timeseries2 = timeseries2.data

    if iscell is not None:
        iscell = iscell.data
        ind  = np.where(iscell > 0)[0]
        timeseries1 = timeseries1[ind, :]

    # data shold be time x component matrix
    neural_data = timeseries1['timeseries'].transpose()  # neural data

    behavioral_data = timeseries2[params['behavior']].values.reshape(-1, 1)    # behavioral dataframe

    # preprocessing
    # sc = StandardScaler(
    #     with_mean=params['standardization_n_mean'],
    #     with_std=params['standardization_n_std']
    # )
    # neural_data = sc.fit_transform(neural_data)
    from utils import standard_norm
    neural_data = standard_norm(neural_data, params['standard_n_mean'], params['standard_n_std'])

    # sc = StandardScaler(
    #     with_mean=params['standardization_b_mean'],
    #     with_std=params['standardization_b_std']
    # )
    # behavioral_data = sc.fit_transform(behavioral_data)
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
    info['actual_predicted'] = TimeSeriesData(
        np.array([Res._endog, Res.mu]).transpose(),
        func_name='glm',
        file_name='actual_predicted'
    )

    return info
