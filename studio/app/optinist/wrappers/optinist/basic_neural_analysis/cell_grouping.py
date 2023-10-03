from studio.app.common.dataclass import TimeSeriesData
from studio.app.optinist.dataclass import NWBFile


def cell_grouping(
    neural_data: TimeSeriesData,
    output_dir: str,
    nwbfile: NWBFile = None,
    params: dict = None,
    **kwargs,
) -> dict():
    import numpy as np

    neural_data = neural_data.data
    std = neural_data.std

    if params["transpose"]:
        neural_data = neural_data.transpose()

    baseline = np.mean(neural_data, axis=1, keepdims=True)
    std = neural_data[np.argmax(neural_data, axis=1)].std()

    grouped_cells = (neural_data.max(axis=1) - baseline) / std > params["threshold"]

    info = {}
    info["grouped_cells"] = TimeSeriesData(
        neural_data[grouped_cells], std=std, file_name="grouped_cells_mean"
    )

    return info
