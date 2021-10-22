def suite2p_roi(ops, opts=None):
	import numpy as np
	from suite2p import extraction, classification, detection, ROI

	# ROI detection
	classfile = classification.user_classfile
	ops, stat = detection.detect(ops=ops, classfile=classfile)

	######## ROI EXTRACTION ##############
	ops, stat, F, Fneu, F_chan2, Fneu_chan2 = extraction.create_masks_and_extract(ops, stat)

	######## ROI CLASSIFICATION ##############
	iscell = classification.classify(stat=stat, classfile=classfile)
	iscell = iscell[:, 0].astype(bool)

	arrays = []
	for i, s in enumerate(stat):
		array = ROI(
			ypix=s['ypix'], xpix=s['xpix'], lam=s['lam'], med=s['med'], do_crop=False
		).to_array(Ly=ops['Ly'], Lx=ops['Lx'])
		array *= i + 1
		arrays.append(array)

	im = np.stack(arrays)
	im[im == 0] = np.nan

	ops['ROI_found'] = np.nanmax(im, axis=0)
	ops['non_cell_roi'] = np.nanmax(im[~iscell])
	ops['cell_roi'] = np.nanmax(im[iscell], axis=0)
	ops['F'] = F
	ops['Fneu'] = Fneu

	ops['max_proj'] = ops['max_proj']
	ops['images'] = ops['max_proj']

	return ops
