#  decoding neural activity by GLM
#  input:  A:matrix[num_cell x timeseries ]   B:timeseries(behavior)[1 x timeseries]
#  Generalized linear modeling using statsmodels
#
#　https://www.statsmodels.org/stable/glm.html

from wrappers.data_wrapper import *
from wrappers.args_check import args_check

@args_check
def GLM(timeseries1: TimeSeriesData, timeseries2: TimeSeriesData,  iscell: IscellData, params: dict):

    # modules specific to function
    import statsmodels.api as sm
    from sklearn.preprocessing import StandardScaler
    import pandas as pd

    timeseries1 = timeseries1.data
    iscell = iscell.data
    ind  = np.where(iscell > 0)[0]
    timeseries1 = timeseries1[ind, :]

    # data shold be time x component matrix
    N = timeseries1['timeseries'].transpose()  # neural data

    B = timeseries2[params['behavior']].values.reshape(-1, 1)    # behavioral dataframe

    # preprocessing  ##################
    sc = StandardScaler(with_mean = params['standardization_n_mean'], with_std = params['standardization_n_std'])
    tN = sc.fit_transform(N)

    sc = StandardScaler(with_mean = params['standardization_b_mean'], with_std = params['standardization_b_std'])
    tB = sc.fit_transform(B)

    # calculate GLM ##################
    tN = pd.DataFrame(tN)
    tB = pd.DataFrame(tB)

    if(params['add_constant']):
        tN = sm.add_constant(tN, prepend=False)

    # set link function
    link = getattr(sm.genmod.families.links, params['link'])()

    # set family
    family = eval('sm.families.' + params['family'] + '(link=link)')

    # model fit
    model = sm.GLM(
        tB, tN, family=family, offset=params['offset'], 
        exposure=params['exposure'], missing=params['missing'])
    Res = model.fit()

    # output structures ###################
    Out_nwb = {}

    # main results
    Out_nwb['params'] = Res.params  # fitted coefficients
    Out_nwb['pvalues'] = Res.pvalues
    Out_nwb['tvalues'] = Res.tvalues # z
    # Out_nwb['stderr'] =

    Out_nwb['aic'] = Res.aic
    Out_nwb['bic'] = Res.bic
    Out_nwb['llf'] = Res.llf # log-Likelihood
    Out_nwb['pearson_chi2'] = Res.pearson_chi2
    Out_nwb['df_model'] = Res.df_model
    Out_nwb['df_resid'] = Res.df_resid
    Out_nwb['scale'] = Res.scale
    Out_nwb['mu'] = Res.mu
    Out_nwb['endog'] = Res._endog

    # additional results
    Out_additional = {'GLM':Res}

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
