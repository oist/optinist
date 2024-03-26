import numpy as np

from studio.app.common.dataclass import ImageData
from studio.app.optinist.core.nwb.nwb import NWBDATASET
from studio.app.optinist.dataclass import EditRoiData, FluoData, IscellData, RoiData


def lccd_detect(
    mc_images: ImageData, output_dir: str, params: dict = None, **kwargs
) -> dict(fluorescence=FluoData, cell_roi=RoiData):
    from studio.app.optinist.wrappers.lccd.lccd_python.lccd import LCCD

    function_id = output_dir.split("/")[-1]
    print("start lccd_detect:", function_id)

    print("params: ", params)
    lccd = LCCD(params)
    D = LoadData(mc_images)
    assert len(D.shape) == 3, "input array should have dimensions (width, height, time)"
    roi = lccd.apply(D)

    dff_f0_frames = params["dff"]["f0_frames"]
    dff_f0_percentile = params["dff"]["f0_percentile"]
    num_cell = roi.shape[1]
    num_frames = D.shape[2]
    iscell = np.ones(num_cell, dtype=int)

    reshapedD = D.reshape([D.shape[0] * D.shape[1], D.shape[2]])
    timeseries = np.zeros([num_cell, num_frames])
    roi_list = []

    for i in range(num_cell):
        roi_list.append((roi[:, i].reshape([D.shape[0], D.shape[1]])) * (i + 1))
        timeseries[i, :] = np.mean(reshapedD[roi[:, i] > 0, :], axis=0)

    im = np.stack(roi_list)
    im[im == 0] = np.nan
    im -= 1

    timeseries_dff = np.ones([num_cell, num_frames]) * np.nan
    for i in range(num_cell):
        for k in range(num_frames):
            if (k - dff_f0_frames >= 0) and (k + dff_f0_frames < num_frames):
                f0 = np.percentile(
                    timeseries[i, k - dff_f0_frames : k + dff_f0_frames],
                    dff_f0_percentile,
                )
                timeseries_dff[i, k] = (timeseries[i, k] - f0) / f0

    roi_list = [{"image_mask": roi[:, i].reshape(D.shape[:2])} for i in range(num_cell)]

    nwbfile = {}
    nwbfile[NWBDATASET.ROI] = {function_id: roi_list}
    nwbfile[NWBDATASET.POSTPROCESS] = {function_id: {"all_roi_img": im}}

    nwbfile[NWBDATASET.COLUMN] = {
        function_id: {
            "name": "iscell",
            "description": "two columns - iscell & probcell",
            "data": iscell,
        }
    }

    nwbfile[NWBDATASET.FLUORESCENCE] = {
        function_id: {
            "Fluorescence": {
                "table_name": "Fluorescence",
                "region": list(range(len(timeseries))),
                "name": "Fluorescence",
                "data": timeseries,
                "unit": "lumens",
            }
        }
    }

    info = {
        "cell_roi": RoiData(
            np.nanmax(im[iscell != 0], axis=0),
            output_dir=output_dir,
            file_name="cell_roi",
        ),
        "fluorescence": FluoData(timeseries, file_name="fluorescence"),
        "dff": FluoData(timeseries_dff, file_name="dff"),
        "iscell": IscellData(iscell),
        "edit_roi_data": EditRoiData(images=mc_images.data, im=im),
        "nwbfile": nwbfile,
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
