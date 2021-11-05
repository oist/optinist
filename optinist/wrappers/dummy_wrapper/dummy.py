import numpy as np
import sys

from wrappers.data_wrapper import *
from wrappers.args_check import args_check

@args_check
def dummy_image2image(image: ImageData, params: dict=None):
    
    """
        get image
        return image
    """
    info = {}
    info['image2image'] = ImageData(np.random.rand((100)).reshape(10, 10), sys._getframe().f_code.co_name)
    return info


@args_check
def dummy_image2time(image: ImageData, params: dict=None):
    
    """
        get image
        return image
    """
    info = {}
    info['image2time'] = TimeSeriesData(np.random.rand((10000)).reshape(10, 1000), sys._getframe().f_code.co_name)
    return info


@args_check
def dummy_image2heat(image: ImageData, params: dict=None):
    
    """
        get image
        return image
    """
    info = {}
    info['image2heat'] = CorrelationData(np.random.rand((10000)).reshape(100, 100), sys._getframe().f_code.co_name)
    return info


@args_check
def dummy_time2time(timeseries: TimeSeriesData, params: dict=None):
    
    """
        get image
        return image
    """
    info = {}
    info['time2time'] = TimeSeriesData(np.random.rand((10000)).reshape(10, 1000), sys._getframe().f_code.co_name)
    return info


@args_check
def dummy_image2image8time(image1: ImageData, params: dict=None):
    
    """
        get image
        return image
    """
    info = {}
    info['image'] = ImageData(np.random.rand((100)).reshape(10, 10), sys._getframe().f_code.co_name)
    info['timeseries'] = TimeSeriesData(np.random.rand((10000)).reshape(10, 1000), sys._getframe().f_code.co_name)
    return info


@args_check
def dummy_image8image2image8time(image1: ImageData, image2: ImageData, params: dict=None):
    
    """
        get image
        return image
    """
    info = {}
    info['image'] = ImageData(np.random.rand((100)).reshape(10, 10), sys._getframe().f_code.co_name)
    info['timeseries'] = TimeSeriesData(np.random.rand((10000)).reshape(10, 1000), sys._getframe().f_code.co_name)
    return info


@args_check
def dummy_time8image2image8time(timeseries: TimeSeriesData, image: ImageData, params: dict=None):
    
    """
        get image
        return image
    """
    info = {}
    info['image'] = ImageData(np.random.rand((100)).reshape(10, 10), sys._getframe().f_code.co_name)
    info['timeseries'] = TimeSeriesData(np.random.rand((100)).reshape(10, 10), sys._getframe().f_code.co_name)
    return info
