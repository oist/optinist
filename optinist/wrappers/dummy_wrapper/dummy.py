import numpy as np
import sys

from wrappers.data_wrapper import *
from wrappers.args_check import args_check
from pynwb import NWBFile

@args_check
def dummy_image2image(
    image: ImageData, nwbfile: NWBFile=None, params: dict=None
    ) -> {'image2image': ImageData}:

    """
        get image
        return image
    """
    info = {}
    info['image2image'] = ImageData(
        np.random.rand((100)).reshape(10, 10),
        func_name=sys._getframe().f_code.co_name)
    return info


@args_check
def dummy_image2time(
    image: ImageData, nwbfile: NWBFile=None, params: dict=None
    ) -> {'image2time': TimeSeriesData}:

    """
        get image
        return image
    """
    info = {}
    info['image2time'] = TimeSeriesData(
        np.random.rand((10000)).reshape(10, 1000),
        func_name=sys._getframe().f_code.co_name)
    return info


@args_check
def dummy_image2heat(
    image: ImageData, nwbfile: NWBFile=None, params: dict=None
    ) -> {'image2heat': CorrelationData}:

    """
        get image
        return image
    """
    info = {}
    info['image2heat'] = CorrelationData(
        np.random.rand((10000)).reshape(100, 100),
        func_name=sys._getframe().f_code.co_name)
    return info


@args_check
def dummy_time2time(
    timeseries: TimeSeriesData, nwbfile: NWBFile=None, params: dict=None
    ) -> {'time2time': TimeSeriesData}:

    """
        get image
        return image
    """
    info = {}
    info['time2time'] = TimeSeriesData(
        np.random.rand((10000)).reshape(10, 1000),
        func_name=sys._getframe().f_code.co_name)
    return info


@args_check
def dummy_image2image8time(
    image1: ImageData, nwbfile: NWBFile=None, params: dict=None
    ) -> {'image': ImageData, 'timeseries': TimeSeriesData}:

    """
        get image
        return image
    """
    info = {}
    info['image'] = ImageData(
        np.random.rand((100)).reshape(10, 10),
        func_name=sys._getframe().f_code.co_name)
    info['timeseries'] = TimeSeriesData(
        np.random.rand((10000)).reshape(10, 1000),
        func_name=sys._getframe().f_code.co_name)
    return info


@args_check
def dummy_image8image2image8time(
    image1: ImageData, image2: ImageData, nwbfile: NWBFile=None, params: dict=None
    ) -> {'image': ImageData, 'timeseries': TimeSeriesData}:

    """
        get image
        return image
    """
    info = {}
    info['image'] = ImageData(
        np.random.rand((100)).reshape(10, 10),
        func_name=sys._getframe().f_code.co_name,
        file_name='image')
    info['timeseries'] = TimeSeriesData(
        np.random.rand((10000)).reshape(10, 1000),
        func_name=sys._getframe().f_code.co_name,
        file_name='timeseries')
    return info


@args_check
def dummy_time8image2image8time(
    timeseries: TimeSeriesData, image: ImageData, nwbfile: NWBFile=None, params: dict=None):

    """
        get image
        return image
    """
    info = {}
    info['image'] = ImageData(
        np.random.rand((100)).reshape(10, 10),
        func_name=sys._getframe().f_code.co_name,
        file_name='image')
    info['timeseries'] = TimeSeriesData(
        np.random.rand((100)).reshape(10, 10),
        func_name=sys._getframe().f_code.co_name,
        file_name='timeseries')
    return info


@args_check
def dummy_keyerror(image: ImageData, nwbfile: NWBFile=None, params: dict=None):

    """
        get image
        return image
    """
    print(params['AAA'])
    info = {}
    return info


@args_check
def dummy_typeerror(image: str, nwbfile: NWBFile=None, params: dict=None):

    """
        get image
        return image
    """
    info = {}
    return info


@args_check
def dummy_image2time8iscell(
    image1: ImageData, nwbfile: NWBFile=None,  params: dict=None
    ) -> {'timeseries': TimeSeriesData, 'iscell': IscellData}:

    """
        get image
        return image
    """
    info = {}
    info['timeseries'] = TimeSeriesData(
        np.random.rand((100)).reshape(10, 10),
        func_name=sys._getframe().f_code.co_name,
        file_name='image')
    info['iscell'] = IscellData(
        np.random.rand((100)).reshape(10, 10),
        func_name=sys._getframe().f_code.co_name,
        file_name='timeseries')
    return info
