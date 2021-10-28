import numpy as np
import sys

from wrappers.data_wrapper import ImageData, TimeSeriesData
from wrappers.args_check import args_check

@args_check
def dummy_image2image(image: ImageData, opts: dict=None):
	
	"""
		get image
		return image
	"""
	info = {}
	info['image2image'] = ImageData(np.random.rand((100)).reshape(10, 10), sys._getframe().f_code.co_name)
	return info


@args_check
def dummy_image2time(image: ImageData, opts: dict=None):
	
	"""
		get image
		return image
	"""
	info = {}
	info['image2time'] = TimeSeriesData(np.random.rand((100)).reshape(10, 10), sys._getframe().f_code.co_name)
	return info


@args_check
def dummy_time2time(timeseries: TimeSeriesData, opts: dict=None):
	
	"""
		get image
		return image
	"""
	info = {}
	info['time2time'] = TimeSeriesData(np.random.rand((100)).reshape(10, 10), sys._getframe().f_code.co_name)
	return info


@args_check
def dummy_image8image2image8time(image1: ImageData, image2: ImageData, opts: dict=None):
	
	"""
		get image
		return image
	"""
	info = {}
	info['image'] = ImageData(np.empty((3, 64, 64)), sys._getframe().f_code.co_name)
	info['timeseries'] = TimeSeriesData(np.empty((100, 100)), sys._getframe().f_code.co_name)
	return info

@args_check
def dummy_time8image2image8time(timeseries: TimeSeriesData, image: ImageData, opts: dict=None):
	
	"""
		get image
		return image
	"""
	info = {}
	info['image'] = ImageData(np.empty((3, 64, 64)), sys._getframe().f_code.co_name)
	info['timeseries'] = TimeSeriesData(np.empty((100, 100)), sys._getframe().f_code.co_name)
	return info
