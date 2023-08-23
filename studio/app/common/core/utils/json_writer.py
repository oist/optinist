import json
import os
from dataclasses import asdict
from typing import Optional

import numpy as np
import pandas as pd
import tifffile

from studio.app.common.core.utils.filepath_creater import (
    create_directory,
    join_filepath,
)
from studio.app.common.schemas.outputs import PlotMetaData


class JsonWriter:
    @classmethod
    def write(cls, filepath, data, orient=None, meta: Optional[PlotMetaData] = None):
        pd.DataFrame(data).to_json(indent=4)
        df = json.loads(pd.DataFrame(data).to_json(indent=4, orient=orient))
        if meta:
            df.update({"meta": asdict(meta)})
        with open(filepath, "w") as f:
            json.dump(df, f, indent=4)

    @classmethod
    def write_as_split(cls, filepath, data, meta: Optional[PlotMetaData] = None):
        cls.write(filepath, data, orient="split", meta=meta)


def save_tiff2json(tiff_filepath, save_dirpath, start_index=None, end_index=None):
    # Tiff画像を読み込む
    tiffs = []
    image = tifffile.imread(tiff_filepath)
    if image.ndim == 2:
        image = image[np.newaxis, :, :]

    for i, page in enumerate(image):
        if i < start_index - 1:
            continue

        if i >= end_index:
            break

        tiffs.append(page.tolist())

    filename, _ = os.path.splitext(os.path.basename(tiff_filepath))
    create_directory(save_dirpath)

    JsonWriter.write_as_split(
        join_filepath(
            [save_dirpath, f"{filename}_{str(start_index)}_{str(end_index)}.json"]
        ),
        pd.DataFrame(tiffs),
    )
