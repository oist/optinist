from optinist.api.dataclass.dataclass import *
from optinist.wrappers.optinist_wrapper.utils import standard_norm
from optinist.api.nwb.nwb import NWBDATASET

def TSNE(
        neural_data: FluoData,
        iscell: IscellData=None,
        params: dict=None
    ) -> dict():

    from sklearn.manifold import TSNE

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

    # preprocessing  ##################
    tX = standard_norm(X, params['standard_mean'], params['standard_std'])

    # calculate TSNE ##################
    tsne = TSNE(**params['TSNE'])

    proj_X = tsne.fit_transform(tX)

    # NWB追加
    nwbfile = {}
    nwbfile[NWBDATASET.POSTPROCESS] = {
        'projectedNd': proj_X
    }

    info = {
        'projectedNd': ScatterData(proj_X, file_name='projectedNd'),
        'nwbfile': nwbfile,
    }

    return info
