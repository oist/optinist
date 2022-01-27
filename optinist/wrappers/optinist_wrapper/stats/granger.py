from wrappers.data_wrapper import *
from wrappers.args_check import args_check
from wrappers.optinist_wrapper.utils import standard_norm

@args_check
def Granger(
        timeseries: TimeSeriesData, iscell: IscellData=None, params: dict=None
    ) -> {'Granger_fval_mat': CorrelationData}:

    # modules specific to function
    # from sklearn.preprocessing import StandardScaler
    from statsmodels.tsa.stattools import grangercausalitytests, adfuller, coint
    import itertools
    from tqdm import tqdm

    timeseries = timeseries.data

    # data shold be time x component matrix
    if iscell is not None:
        iscell = iscell.data
        ind  = np.where(iscell > 0)[0]
        timeseries = timeseries[ind, :]

    if params['transpose']:
        X = timeseries.transpose()
    else:
        X = timesereies

    num_cell = X.shape[0]
    comb = list(itertools.permutations(range(num_cell), 2))  # combinations with dup
    num_comb= len(comb)

    # # preprocessing  ##################
    tX = standard_norm(X, params['standard_mean'], params['standard_std'])

    # calculate dickey-fuller test  ##################
    # augmented dickey-fuller test---- if p val is large, it cannot reject  there is a unit root
    # -> small p-val means OK : means this is not unit root process that it can apply Causality test

    adf = {}
    adf['adf_teststat'] = np.zeros([num_cell], dtype = 'float64')
    adf['adf_pvalue'] = np.zeros([num_cell], dtype = 'float64')
    adf['adf_usedlag'] = np.zeros([num_cell], dtype = 'int')
    adf['adf_nobs'] = np.zeros([num_cell], dtype = 'int')
    adf['adf_critical_values'] = np.zeros([num_cell, 3], dtype = 'float64')
    adf['adf_icbest'] = np.zeros([num_cell], dtype = 'float64')

    if (params['adfuller_test']):
        print('adfuller test ')

        for i in tqdm(range(num_cell)):
            tp = adfuller(
                tX[:, i],
                maxlag=params['adfuller_maxlag'],
                regression=params['adfuller_regression'],
                autolag=params['adfuller_autolag'],
                store=params['adfuller_store'],
                regresults=params['adfuller_regresults']
            )

            adf['adf_teststat'][i] = tp[0]
            adf['adf_pvalue'][i] = tp[1]
            adf['adf_usedlag'][i] = tp[2]
            adf['adf_nobs'][i] = tp[3]
            adf['adf_critical_values'][i, :] = np.array([tp[4]['1%'],  tp[4]['5%'], tp[4]['10%']])
            adf['adf_icbest'][i] = tp[5]

    #  test for cointegration ##################
    # augmented engle-granger two-step test ----- Test for no-cointegration of a univariate equation
    # if p val is small, the relation is cointegration
    # -> check this if ADF pval is large
    cit={}
    cit['cit_count_t'] = np.zeros([num_comb], dtype='float64')
    cit['cit_pvalue'] = np.zeros([num_comb], dtype='int')
    cit['cit_crit_value'] = np.zeros([num_comb, 3], dtype='float64')

    if (params['coint_test']):
        print('cointegration test ')

        for i in tqdm(range(num_comb)):
            tp = coint(
                X[:, comb[i][0]],
                X[:, comb[i][1]],
                trend=params['coint_trend'],
                method=params['coint_method'],
                maxlag=params['coint_maxlag'],
                autolag=params['coint_autolag']
            )
            if not np.isnan(tp[0]):
                cit['cit_count_t'][i] = tp[0]
            if not np.isnan(tp[1]):
                cit['cit_pvalue'][i] = tp[1]
            cit['cit_crit_value'][i, :] = tp[2]

    #  Granger causality ##################
    print('granger test ')

    if( hasattr(params['Granger_maxlag'], '__iter__')):
        num_lag = len(params['Granger_maxlag'])
    else:
        num_lag = params['Granger_maxlag']

    GC = {}
    GC['gc_combinations'] = comb
    GC['gc_ssr_ftest'] = np.zeros([num_comb, num_lag, 4], dtype='float64') # test statistic, pvalues, deg-free
    GC['gc_ssr_chi2test'] = np.zeros([num_comb,  num_lag, 3], dtype='float64')
    GC['gc_lrtest'] = np.zeros([num_comb, num_lag, 3], dtype='float64')
    GC['gc_params_ftest'] = np.zeros([num_comb, num_lag, 4], dtype='float64')
    GC['gc_OLS_restricted'] = [[0] * num_lag for i in range(num_comb)]
    GC['gc_OLS_unrestricted'] = [[0] * num_lag for i in range(num_comb)]
    GC['gc_OLS_restriction_matrix'] = [[0] * num_lag for i in range(num_comb)]

    # for figure presentation
    GC['Granger_fval_mat'] = [np.zeros([num_cell, num_cell]) for i in range(num_lag)]

    for i in tqdm(range(len(comb))):
        # The Null hypothesis for grangercausalitytests is that the time series in the second column1,
        # does NOT Granger cause the time series in the first column0
        # column 1 -> colum 0
        tp = grangercausalitytests(
            X[:, [comb[i][0], comb[i][1]]],
            params['Granger_maxlag'],
            verbose=False,
            addconst=params['Granger_addconst']
        )

        for j in range(len(tp)): # number of lag
            GC['gc_ssr_ftest'][i, j, :] = tp[j+1][0]['ssr_ftest'][0:4]  #ssr based F test (F, pval, df_denom, df_num)
            GC['gc_ssr_chi2test'][i, j, :] = tp[j+1][0]['ssr_chi2test'][0:3] #ssr based chi2test (chi2, pval, df)
            GC['gc_lrtest'][i, j, :] = tp[j+1][0]['lrtest'][0:3] #likelihood ratio test (chi2, pval, df)
            GC['gc_params_ftest'][i, j, :] = tp[j+1][0]['params_ftest'][0:4] #parameter F test (F, pval, df_denom, df_num)
            GC['gc_OLS_restricted'][i][j] = tp[j+1][1][0]
            GC['gc_OLS_unrestricted'][i][j] = tp[j+1][1][1]
            GC['gc_OLS_restriction_matrix'][i][j] = tp[j+1][1][2]

            GC['Granger_fval_mat'][j][comb[i][0], comb[i][1]] = tp[j+1][0]['ssr_ftest'][0]

    GC['Granger_fval_mat'] = np.array(GC['Granger_fval_mat'])

    # output structures ###################
    Out_nwb = {}

    # main results
    Out_nwb['gc_combinations'] = GC['gc_combinations']
    Out_nwb['gc_ssr_ftest'] = GC['gc_ssr_ftest']
    Out_nwb['gc_ssr_chi2test'] = GC['gc_ssr_chi2test']
    Out_nwb['gc_lrtest'] = GC['gc_lrtest']
    Out_nwb['gc_params_ftest'] = GC['gc_params_ftest']

    Out_nwb['cit_pvalue'] = cit['cit_pvalue']
    Out_nwb['adf_pvalue'] = adf['adf_pvalue']

    # additional results
    Out_additional = {}
    Out_additional.update(cit)
    Out_additional.update(adf)
    Out_additional.update(GC)

    # main results for plot
    info = {}
    info['Granger_fval_mat'] = CorrelationData(
        GC['Granger_fval_mat'][0],
        func_name='granger',
        file_name='gfm'
    )

    return info
