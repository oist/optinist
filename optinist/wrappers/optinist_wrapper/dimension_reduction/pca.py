from optinist.wrappers.data_wrapper import *
from optinist.wrappers.optinist_wrapper.utils import standard_norm
from optinist.api.nwb.nwb import NWBDATASET

def PCA(
        neural_data: FluoData,
        iscell: IscellData=None,
        nwbfile: NWBFile=None,
        params: dict=None
    ) -> dict():

    # modules specific to function
    from sklearn.decomposition import PCA

    neural_data = neural_data.data

    # data shold be time x component matrix
    if params['transpose']:
        X = neural_data.transpose()
    else:
        X = neural_data

    if iscell is not None:
        iscell = iscell.data
        ind  = np.where(iscell > 0)[0]
        X = X[:, ind]

    # # preprocessing  ##################
    tX = standard_norm(X, params['standard_mean'], params['standard_std'])

    # calculate PCA ##################
    pca = PCA(**params['PCA'])
    proj_X = pca.fit_transform(tX)

    # NWB追加
    if nwbfile is not None:
        nwbfile[NWBDATASET.POSTPROCESS] = {
            'pca_projectedNd': proj_X,
            'explained_variance': pca.explained_variance_ratio_,
        }

    info = {
        'explained_variance': BarData(pca.explained_variance_ratio_, file_name='evr'),
        'projectedNd': ScatterData(proj_X, file_name='projectedNd'),
        'nwbfile': nwbfile,
    }

    return info
