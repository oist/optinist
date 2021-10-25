from wrappers.data_wrapper import ImageData, TimeSeriesData
import numpy as np
from wrappers.args_check import args_check

@args_check
def dummy_image2image(image: ImageData, opts=None):
	
	"""
		get image
		return image
	"""
	info = {}
	info['timeseries'] = TimeSeriesData(np.empty((3, 64, 64)))
	info['images'] = ImageData(np.empty((3, 64, 64)))
	return info


@args_check
def dummy_image2time(image: ImageData, opts=None):
	
	"""
		get image
		return image
	"""
	info = {}
	info['timeseries'] = TimeSeriesData(np.empty((100, 100)))
	return info


@args_check
def dummy_time2time(timeseries: TimeSeriesData, opts=None):
	
	"""
		get image
		return image
	"""
	info = {}
	info['timeseries'] = TimeSeriesData(np.empty((100, 100)))
	return info

@args_check
def dummy_image8image2image(image1: ImageData, image2: ImageData, opts=None):
	
	"""
		get image
		return image
	"""
	info = {}
	info['timeseries'] = TimeSeriesData(np.empty((100, 100)))
	return info