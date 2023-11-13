import numpy as np
import scipy.sparse

from studio.app.optinist.wrappers.lccd.lccd_python import utils


def filter_roi_by_area(roi, min_area, max_area):
    """
    Parameters
    ----------
    roi: np.array or scipy.sparse.csc (Maybe work with other sparse matrix types)
        roi has shape (X, n).
        X is spatial dimension.
        n is the number of rois
    min_area: int
    max_area: int

    Returns
    -------
    roi
    """

    n_rois = roi.shape[1]
    indices = []
    for i in range(n_rois):
        area = np.sum(roi[:, i])
        if min_area <= area and area <= max_area:
            indices.append(i)
    return roi[:, np.array(indices)]


class RoiIntegration:
    def __init__(self, overlap_threshold, min_area, max_area, sparse=True, **kwargs):
        self.overlap_threshold = overlap_threshold
        self.min_area = min_area
        self.max_area = max_area
        self.sparse = sparse

    def delete_col(self, arr, idx):
        if self.sparse:
            return utils.delete_col_csc(arr, idx)
        else:
            return np.delete(arr, idx, 1)

    def cat_rois(self, roi_a, roi_b):
        if self.sparse:
            return scipy.sparse.hstack((roi_a, roi_b))
        return np.hstack((roi_a, roi_b))

    def get_area(self, roi, i):
        """
        get area of ith region in roi.
        i.e. get sum of ith column of 2d array.
        """
        if self.sparse:
            return utils.count_nonzero_values_in_col_csc(roi, i)
        return np.sum(
            roi[:, i]
        )  # for dense array, np.count_nonzero is slower than np.sum(col)

    def gather_overlapping_regions(self, sim):
        """
        Gather regions those have overlapping area with other regions.
        """
        if self.sparse:
            return [
                region_in_a
                for region_in_a in range(sim.shape[0])
                if sim.getrow(region_in_a).nnz > 0
            ]
        return [
            region_in_a
            for region_in_a in range(sim.shape[0])
            if np.sum(sim[region_in_a]) > 0
        ]

    def remove_overlap(self, roi_a, roi_b, i, j):
        """
        This is inplace operation.
        From both roi_a's ith region and roi_b's jth region,
        remove their intersection.

        Let
        roi_a[:, i] = [1, 0, 1, 1, 0, 1, 1].T
        roi_b[:, j] = [0, 1, 1, 0, 1, 1, 0].T

        intersection is [0, 0, 1, 0, 0, 1, 0].T
        then, calling remove_overlap(roi_a, roi_b, i, j),
        roi_a[:, i] = [1, 0, 0, 1, 0, 0, 1].T
        roi_b[:, j] = [0, 1, 0, 0, 1, 0, 0].T
        """

        if self.sparse:  # TODO: improve readability
            roi_a_region_indices = roi_a.indices[roi_a.indptr[i] : roi_a.indptr[i + 1]]
            roi_b_region_indices = roi_b.indices[roi_b.indptr[j] : roi_b.indptr[j + 1]]
            overlap_indices = utils.intersect_unique_sorted_1d(
                roi_a_region_indices, roi_b_region_indices
            )
            roi_a_region_separated_indices = np.setdiff1d(
                roi_a_region_indices, overlap_indices, True
            )  # This set operation may be slow.
            roi_b_region_separated_indices = np.setdiff1d(
                roi_b_region_indices, overlap_indices, True
            )
            utils.delete_col_csc_inplace(roi_a, i)
            utils.insert_binary_col_binary_csc(roi_a, roi_a_region_separated_indices, i)
            utils.delete_col_csc_inplace(roi_b, j)
            utils.insert_binary_col_binary_csc(roi_b, roi_b_region_separated_indices, j)
        else:
            region_non_overlap = (roi_a[:, i] * roi_b[:, j]) == 0
            region_non_overlap = region_non_overlap.astype(roi_a.dtype)
            roi_a[:, i] *= region_non_overlap
            roi_b[:, j] *= region_non_overlap

    def apply(self, roi_a, roi_b):
        # requires test
        if roi_b is None:
            return roi_a
        dtype = np.float32
        roi_a = roi_a.astype(dtype)
        roi_b = roi_b.astype(dtype)
        if self.sparse:
            roi_a = scipy.sparse.csc_matrix(roi_a)
            roi_b = scipy.sparse.csc_matrix(roi_b)
        while 1:
            change_occured = False
            sim = roi_a.T.dot(roi_b)
            # sim[i][j] is area of overlap
            # by ith region in roi_a and jth region in roi_b
            # roi_a.shape = (1500x1500, 5000), roi_b.shape = (1500x1500, 3000)
            # density 30 / (1500x1500)
            # roi_a.T.dot(roi_b) took 6ms (AMD Ryzen 5900x)

            overlapping_regions_in_a = self.gather_overlapping_regions(sim)
            for region_in_roi_a in overlapping_regions_in_a:
                region_in_roi_b = np.argmax(
                    sim[region_in_roi_a]
                )  # most overlapping region in roi_b
                area_overlap = sim[region_in_roi_a, region_in_roi_b]
                if area_overlap == 0:
                    continue
                area_a = self.get_area(roi_a, region_in_roi_a)
                area_b = self.get_area(roi_b, region_in_roi_b)
                if (
                    area_overlap / area_a > self.overlap_threshold
                    or area_overlap / area_b > self.overlap_threshold
                ):
                    change_occured = True
                    if area_a > area_b:
                        roi_b = self.delete_col(roi_b, region_in_roi_b)
                        # column deletion from large array is slow.
                    else:
                        roi_a = self.delete_col(roi_a, region_in_roi_a)
                    break
                elif (
                    area_overlap / area_a <= self.overlap_threshold
                    and area_overlap / area_b <= self.overlap_threshold
                ):
                    change_occured = True
                    self.remove_overlap(roi_a, roi_b, region_in_roi_a, region_in_roi_b)
                    break
                else:
                    continue
            if not change_occured:
                break
        roi = self.cat_rois(roi_a, roi_b)
        roi = filter_roi_by_area(roi, self.min_area, self.max_area)
        return roi
