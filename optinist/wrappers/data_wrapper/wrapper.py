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
            self.path = os.path.join('files', func_name, 'images')

            if not os.path.exists(self.path):
                os.makedirs(self.path, exist_ok=True)

            if len(self.data.shape) == 2:
                self.data = self.data[np.newaxis, :, :]

            print(self.data.shape)
            for i in range(len(self.data)):
                img = Image.fromarray(np.uint8(self.data[i]))
                img.save(os.path.join(self.path, f'{str(i)}.png'))



class TimeSeriesData:
    def __init__(self, data, func_name='timeseries'):
        self.data = data
        _dir = os.path.join('files', func_name)

        if not os.path.exists(_dir):
            os.makedirs(_dir, exist_ok=True)

        self.path = os.path.join(_dir, 'timeseries.json')

        pd.DataFrame(self.data).to_json(self.path, indent=4)


class Suite2pData:
    def __init__(self, data):
        self.data = data
