import numpy as np

from studio.app.common.dataclass import (
    HeatMapData,
    ImageData,
    ScatterData,
    TimeSeriesData,
)
from studio.app.optinist.dataclass import IscellData, RoiData


def dummy_image2image(
    image: ImageData, params: dict = None
) -> dict(image2image=ImageData):
    """
    get image
    return image
    """
    info = {}
    info["image2image"] = ImageData(np.random.rand((300_0000)).reshape(300, 100, 100))
    return info


def dummy_image2time(
    image: ImageData, params: dict = None
) -> dict(image2time=TimeSeriesData):
    """
    get image
    return image
    """
    info = {}
    info["image2time"] = TimeSeriesData(np.random.rand((100000)).reshape(100, 1000))
    return info


def dummy_image2heat(
    image: ImageData, params: dict = None
) -> dict(image2heat=HeatMapData):
    """
    get image
    return image
    """
    info = {}
    info["image2heat"] = HeatMapData(np.random.rand((10000)).reshape(100, 100))
    return info


def dummy_time2time(
    timeseries: TimeSeriesData, params: dict = None
) -> dict(time2time=TimeSeriesData):
    """
    get image
    return image
    """
    info = {}
    info["time2time"] = TimeSeriesData(np.random.rand((10000)).reshape(10, 1000))
    return info


def dummy_image2image8time(
    image1: ImageData, params: dict = None
) -> dict(image=ImageData, timeseries=TimeSeriesData):
    """
    get image
    return image
    """
    info = {}
    info["image"] = ImageData(np.random.rand((100)).reshape(10, 10))
    info["timeseries"] = TimeSeriesData(np.random.rand((10000)).reshape(10, 1000))
    return info


def dummy_image8image2image8time(
    image1: ImageData, image2: ImageData, params: dict = None
) -> dict(image=ImageData, timeseries=TimeSeriesData):
    """
    get image
    return image
    """
    info = {}
    info["image"] = ImageData(np.random.rand((100)).reshape(10, 10), file_name="image")
    info["timeseries"] = TimeSeriesData(
        np.random.rand((10000)).reshape(10, 1000), file_name="timeseries"
    )
    return info


def dummy_time8image2image8time(
    timeseries: TimeSeriesData, image: ImageData, params: dict = None
) -> dict():
    """
    get image
    return image
    """
    info = {}
    info["image"] = ImageData(np.random.rand((100)).reshape(10, 10), file_name="image")
    info["timeseries"] = TimeSeriesData(
        np.random.rand((100)).reshape(10, 10), file_name="timeseries"
    )
    return info


def dummy_keyerror(image: ImageData, params: dict = None) -> dict():
    """
    get image
    return image
    """
    print(params["AAA"])
    info = {}
    return info


def dummy_typeerror(image: str, params: dict = None) -> dict():
    """
    get image
    return image
    """
    info = {}
    return info


def dummy_image2time8iscell(
    image1: ImageData, params: dict = None
) -> dict(timeseries=TimeSeriesData, iscell=IscellData):
    """
    get image
    return image
    """
    info = {}
    info["timeseries"] = TimeSeriesData(
        np.random.rand((100)).reshape(10, 10), file_name="image"
    )
    info["iscell"] = IscellData(
        np.random.rand((100)).reshape(10, 10), file_name="timeseries"
    )
    return info


def dummy_image2roi(image1: ImageData, params: dict = None) -> dict(roi=RoiData):
    """
    get image
    return image
    """
    roi_data = np.random.randint(low=1, high=1000, size=100_00).astype(float)
    import random

    null_data = np.array(random.sample(list(np.arange(100_00)), 90_00))
    roi_data[null_data] = np.nan
    roi_data = roi_data.reshape(100, 100)
    info = {}
    info["roi"] = RoiData(roi_data)
    return info


def dummy_image2image8roi(
    image1: ImageData, params: dict = None
) -> dict(image=ImageData, roi=RoiData):
    """
    get image
    return image
    """
    import random

    info = {}
    info["images"] = ImageData(np.random.rand((100_00)).reshape(1, 100, 100))
    roi_data = np.random.rand((100_00))
    null_data = np.array(random.sample(list(np.arange(100_00)), 90_00))
    roi_data[null_data] = np.nan
    roi_data = roi_data.reshape(100, 100)
    info["roi"] = RoiData(roi_data)
    return info


def dummy_image2image8roi8time8heat(
    image1: ImageData, params: dict = None
) -> dict(image=ImageData, roi=RoiData, timeseries=TimeSeriesData, heat=HeatMapData):
    """
    get image
    return image
    """
    import random

    info = {}
    info["images"] = ImageData(np.random.rand((100_00)).reshape(1, 100, 100))
    roi_data = np.random.rand((100_00))
    null_data = np.array(random.sample(list(np.arange(100_00)), 90_00))
    roi_data[null_data] = np.nan
    roi_data = roi_data.reshape(100, 100)
    info["roi"] = RoiData(roi_data)
    info["timeseries"] = TimeSeriesData(
        np.random.rand((100)).reshape(10, 10), file_name="image"
    )
    info["heat"] = HeatMapData(np.random.rand((10000)).reshape(100, 100))
    return info


def dummy_image2scatter(
    image: ImageData, params: dict = None
) -> dict(scatter=ScatterData):
    """
    get image
    return scatter
    """
    info = {}
    info["image2scatter"] = ScatterData(np.random.rand((1000)).reshape(500, 2))
    return info
