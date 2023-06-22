import numpy as np

from optinist.api.dataclass.dataclass import (
    BehaviorData,
    FluoData,
    HeatMapData,
    IscellData,
    TimeSeriesData,
)
from optinist.api.nwb.nwb import NWBDATASET


def calc_trigger(behavior_data, trigger_type, trigger_threshold):
    flg = np.array(behavior_data > trigger_threshold, dtype=int)
    if trigger_type == "up":
        trigger_idx = np.where(np.ediff1d(flg) == 1)[0]
    elif trigger_type == "down":
        trigger_idx = np.where(np.ediff1d(flg) == -1)[0]
    elif trigger_type == "cross":
        trigger_idx = np.where(np.ediff1d(flg) != 0)[0]
    else:
        trigger_idx = np.where(np.ediff1d(flg) == 0)[0]

    return trigger_idx


def calc_trigger_average(neural_data, trigger_idx, start_time, end_time):
    num_frame = neural_data.shape[0]

    ind = np.array(range(start_time, end_time), dtype=int)

    event_trigger_data = []
    for trigger in trigger_idx:
        target_idx = ind + trigger

        if np.min(target_idx) >= 0 and np.max(target_idx) < num_frame:
            event_trigger_data.append(neural_data[target_idx])

    # (num_event, cell_number, event_time_lambda)
    event_trigger_data = np.array(event_trigger_data)

    return event_trigger_data


def ETA(
    neural_data: FluoData,
    behaviors_data: BehaviorData,
    output_dir: str,
    iscell: IscellData = None,
    params: dict = None,
) -> dict(mean=TimeSeriesData):
    neural_data = neural_data.data
    behaviors_data = behaviors_data.data

    if params["transpose_x"]:
        X = neural_data.transpose()
    else:
        X = neural_data

    if params["transpose_y"]:
        Y = behaviors_data.transpose()
    else:
        Y = behaviors_data

    assert (
        X.shape[0] == Y.shape[0]
    ), f"""
        neural_data and behaviors_data is not same dimension,
        neural.shape{X.shape}, behavior.shape{Y.shape}"""

    if iscell is not None:
        iscell = iscell.data
        ind = np.where(iscell > 0)[0]
        cell_numbers = ind
        X = X[:, ind]

    Y = Y[:, params["target_index"]]

    # calculate Triggers
    trigger_idx = calc_trigger(Y, params["trigger_type"], params["trigger_threshold"])

    # calculate Triggered average
    event_trigger_data = calc_trigger_average(
        X, trigger_idx, params["start_time"], params["end_time"]
    )

    # (cell_number, event_time_lambda)
    if len(event_trigger_data) > 0:
        mean = np.mean(event_trigger_data, axis=0)
        sem = np.std(event_trigger_data, axis=0) / np.sqrt(len(event_trigger_data))
        mean = mean.transpose()
        sem = sem.transpose()
    else:
        assert False, "Output data size is 0"

    # NWB追加
    nwbfile = {}
    nwbfile[NWBDATASET.POSTPROCESS] = {
        "mean": mean,
        "sem": sem,
        "num_sample": [len(mean)],
    }

    min_value = np.min(mean, axis=1, keepdims=True)
    max_value = np.max(mean, axis=1, keepdims=True)
    norm_mean = (mean - min_value) / (max_value - min_value)

    info = {}
    info["mean"] = TimeSeriesData(
        mean,
        std=sem,
        index=list(np.arange(params["start_time"], params["end_time"])),
        cell_numbers=cell_numbers if iscell is not None else None,
        file_name="mean",
    )
    info["mean_heatmap"] = HeatMapData(
        norm_mean,
        columns=list(np.arange(params["start_time"], params["end_time"])),
        file_name="mean_heatmap",
    )
    info["nwbfile"] = nwbfile

    return info
