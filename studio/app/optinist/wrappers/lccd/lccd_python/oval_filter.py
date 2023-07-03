import numpy as np
import scipy.sparse
import skimage

from studio.app.optinist.wrappers.lccd.lccd_python import utils


def oval_filter(label_mat, sparse=False):
    """ovel_filter functionality.
    Let label_mat has shape(X, Y),
    returned array has shape (X x Y, number of cells)

    when sparse=True, this method took 137ms
    when sparse=False, this method took 325ms,
    sparse=True gave 2.5x speedup.
    spatial size was 506 x 506, number of cells after filtering was 375.
    """
    if sparse:
        rows, data = [], []
    else:
        rois = []
    regionprops = skimage.measure.regionprops(label_mat)
    for regionprop in regionprops:
        area_ratio = (
            regionprop.axis_major_length
            * regionprop.axis_minor_length
            * np.pi
            / 4
            / regionprop.area
        )
        eccen = regionprop.eccentricity
        if area_ratio < 1.8 and eccen < 0.99:  # TODO Make configurable
            roi = (label_mat == regionprop.label).astype(np.uint8)
            if sparse:
                datum, row = utils.array_to_lil_row(roi.ravel())
                data.append(datum)
                rows.append(row)
            else:
                roi = roi.reshape(-1, 1)
                rois.append(roi)

    if sparse:
        if len(data) == 0:
            raise ValueError("No roi region found. Revise config and input.")
        lil_arr = scipy.sparse.lil_matrix(
            (len(data), label_mat.shape[0] * label_mat.shape[1]), dtype=data[0][0].dtype
        )
        lil_arr.data = np.array(data, dtype="object")
        lil_arr.rows = np.array(rows, dtype="object")
        csr_arr = scipy.sparse.csr_matrix(lil_arr)
        return csr_arr.T  # By transpose, csr_matrix is converted to csc_matrix.
    if len(rois) == 0:
        raise ValueError("No roi region found. Revise config and input.")
    rois = np.concatenate(rois, 1)
    return rois
