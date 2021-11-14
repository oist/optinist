import os
import numpy as np
import pandas as pd
import imageio
from PIL import Image

class ImageData:
    def __init__(self, data, func_name='image'):
        if type(data) == str:
            self.path = data
            self.data = np.array(imageio.volread(self.path))
        else:
            self.data = data
            # self.path = os.path.join('files', func_name, 'images')
            _dir = os.path.join('files', func_name)
            self.path = os.path.join(_dir, 'image.json')

            if not os.path.exists(_dir):
                os.makedirs(_dir, exist_ok=True)

            if len(self.data.shape) == 2:
                self.data = self.data[np.newaxis, :, :]

            # self.data = np.uint8(self.data)
            # import pdb; pdb.set_trace()
            images = []
            for i, _img in enumerate(self.data[:10]):
                images.append(_img.tolist())

            pd.DataFrame(images).to_json(self.path, indent=4, orient="values")



class TimeSeriesData:
    def __init__(self, data, func_name='timeseries'):
        self.data = data
        _dir = os.path.join('files', func_name)

        if not os.path.exists(_dir):
            os.makedirs(_dir, exist_ok=True)

        self.path = os.path.join(_dir, 'timeseries.json')

        pd.DataFrame(self.data).to_json(self.path, indent=4)


class CorrelationData:
    def __init__(self, data, func_name='heatmap'):
        self.data = data

        _dir = os.path.join('files', func_name)

        if not os.path.exists(_dir):
            os.makedirs(_dir, exist_ok=True)

        self.path = os.path.join(_dir, 'correlation.json')

        pd.DataFrame(self.data).to_json(self.path, indent=4, orient="values")


class Suite2pData:
    def __init__(self, data):
        self.data = data


class IscellData:
    def __init__(self, data):
        self.data = data


        