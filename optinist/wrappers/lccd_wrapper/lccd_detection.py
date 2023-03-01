from optinist.api.dataclass.dataclass import *
import numpy as np


def lccd_detect(
    mc_images: ImageData, params: dict = None
) -> dict(fluorescence=FluoData, cell_roi=RoiData):
    from .lccd_python.lccd import LCCD

    print('params: ', params)
    lccd = LCCD(params)
    D = LoadData(mc_images)
    assert len(D.shape) == 3, "input array should have dimensions (width, height, time)"
    roi = lccd.apply(D)

    dff_f0_frames = params['dff']['f0_frames']
    dff_f0_percentile = params['dff']['f0_percentile']
    num_cell = roi.shape[1]
    num_frames = D.shape[2]

    reshapedD = D.reshape([D.shape[0] * D.shape[1], D.shape[2]])
    timeseries = np.zeros([num_cell, num_frames])
    roi_list = []

    for i in range(num_cell):
        roi_list.append((roi[:, i].reshape([D.shape[0], D.shape[1]])) * (i + 1))
        timeseries[i, :] = np.mean(reshapedD[roi[:, i] > 0, :], axis=0)

    roi_list = np.stack(roi_list)
    roi_list[roi_list == 0] = np.nan

    timeseries_dff = np.ones([num_cell, num_frames]) * np.nan
    for i in range(num_cell):
        for k in range(num_frames):
            if (k - dff_f0_frames >= 0) and (k + dff_f0_frames < num_frames):
                f0 = np.percentile(
                    timeseries[i, k - dff_f0_frames : k + dff_f0_frames],
                    dff_f0_percentile,
                )
                timeseries_dff[i, k] = (timeseries[i, k] - f0) / f0

    info = {
        'rois': RoiData(np.nanmax(roi_list, axis=0), file_name='cell_roi'),
        'fluorescence': FluoData(timeseries, file_name='fluorescence'),
        'dff': FluoData(timeseries_dff, file_name='dff'),
    }

    return info


def LoadData(mc_images):
    from PIL import Image

    img_pile = Image.open(mc_images.path[0])
    num_page = img_pile.n_frames

    imgs = np.zeros((img_pile.size[0], img_pile.size[1], num_page))
    for i in range(num_page):
        img_pile.seek(i)
        imgs[:, :, i] = np.asarray(img_pile)
    img_pile.close()

    maxval = np.max(imgs)
    minval = np.min(imgs)
    imgs = (imgs - minval) / (maxval - minval)

    return imgs
