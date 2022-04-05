from wrappers.data_wrapper import *
from wrappers.optinist_wrapper.utils import standard_norm
from wrappers.nwb_wrapper.const import NWBDATASET

def TSNE(
        neural_data: FluoData,
        iscell: IscellData=None,
        nwbfile: NWBFile=None,
        params: dict=None
    ) -> {}:

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

    info = {}
    info['projectedNd'] = ScatterData(
        proj_X,
        file_name='projectedNd'
    )

    # NWB追加
    if nwbfile is not None:
        nwbfile[NWBDATASET.POSTPROCESS] = {
            'projectedNd': proj_X
        }

    info['nwbfile'] = nwbfile

    return info
