from optinist.wrappers.data_wrapper import *
from optinist.wrappers.optinist_wrapper.utils import standard_norm
from optinist.wrappers.nwb_wrapper.const import NWBDATASET

def CCA(
        neural_data: FluoData,
        behaviors_data: BehaviorData,
        iscell: IscellData=None,
        nwbfile: NWBFile=None,
        params: dict=None
    ) -> dict():

    from sklearn.cross_decomposition import CCA

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

    Y = Y[:, params['target_index']].reshape(-1, 1)

    # # preprocessing  ##################
    tX = standard_norm(X, params['standard_x_mean'], params['standard_x_std'])
    tY = standard_norm(Y, params['standard_y_mean'], params['standard_y_std'])

    # calculate CCA 
    cca = CCA(**params['CCA'])
    projX, projY = cca.fit_transform(tX, tY)

    proj = np.concatenate([projX, projY], axis=1)

    info = {}
    info['projectedNd'] = ScatterData(
        proj,
        file_name='projectedNd'
    )

    # NWB追加
    if nwbfile is not None:
        nwbfile[NWBDATASET.POSTPROCESS] = {
            'projectedNd': proj
        }

    info['nwbfile'] = nwbfile

    return info
