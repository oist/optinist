# import argparse
# import joblib
# import json

import numpy as np

from optinist.wrappers.lccd_wrapper.lccd_python.blob_detector import BlobDetector
from optinist.wrappers.lccd_wrapper.lccd_python.roi_integration import RoiIntegration


class LCCD:
    def __init__(self, config):
        self.blob_detector = BlobDetector(**config["blob_detector"])
        self.roi_integration = RoiIntegration(**config["roi_integration"])
        self.frame_divider = config["lccd"]["frame_divider"]

    def apply(self, v):
        # lazy loading of v will make this method memory efficient.
        # https://qiita.com/Hanjin_Liu/items/7a01c1c481161fe20a3e
        roi = None
        X, Y, T = v.shape
        for m in range(0, T // self.frame_divider - 1):
            tm = v[:, :, m * self.frame_divider : (m + 1) * self.frame_divider]
            roi_tm = self.blob_detector.apply(tm)
            roi = self.roi_integration.apply(roi_tm, roi)
        if self.roi_integration.sparse and roi is not None:
            roi = np.array(roi.todense())
        return roi


# def main(args):
#     # args = parser.parse_args()
#     with open("config/config.json", "r") as f:
#         config = json.load(f)
#     lccd = LCCD(config)

#     v = np.load("profiles/src.npy")
#     assert len(v.shape) == 3, (
#         "input array should have"
#         "dimensions (width, height, time)"
#     )
#     roi = lccd.apply(v)
#     joblib.dump({"config": config, "roi": roi}, "out.dump")


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser()
#     parser.add_argument(
#         "--conf", default="config/config.json",
#         help="lccd config file"
#     )
#     parser.add_argument("--src", help="path to input file.")
#     parser.add_argument("--dst", help="path to output file.")
#     args = parser.parse_args()
#     main(args)
