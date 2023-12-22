import numpy as np
from scipy.stats import mode

from studio.app.common.dataclass import HeatMapData, TimeSeriesData
from studio.app.optinist.core.nwb.nwb import NWBDATASET
from studio.app.optinist.dataclass import BehaviorData, FluoData, IscellData


def calc_trigger(behavior_data, trigger_type, trigger_threshold):
    behavior_data = np.array(behavior_data, dtype=float)
    flg = np.array(behavior_data > trigger_threshold, dtype=int)
    trigger_idx = []
    trigger_lengths = []

    def find_trigger_length(index, flag_value):
        length = 0
        while index < len(flg) and flg[index] == flag_value:
            index += 1  # from initial trigger find subsequent end of trigger
            length += 1  # record the length of each trigger
        return length

    i = 0
    while i < len(flg):
        if (
            (trigger_type == "up" and flg[i] == 1 and (i == 0 or flg[i - 1] == 0))
            or (trigger_type == "down" and flg[i] == 0 and i > 0 and flg[i - 1] == 1)
            or (
                trigger_type == "cross"
                and (
                    (flg[i] == 1 and (i == 0 or flg[i - 1] == 0))
                    or (flg[i] == 0 and i > 0 and flg[i - 1] == 1)
                )
            )
        ):
            length = find_trigger_length(i, flg[i])
            trigger_idx.append(i)
            trigger_lengths.append(length)
            i += length
        else:
            i += 1
    if trigger_lengths:
        trigger_len = mode(trigger_lengths).mode[0]  # find most common length
        trigger_idx = [
            idx
            for idx, length in zip(trigger_idx, trigger_lengths)
            if length == trigger_len  # if trigger_len is different (boundary cut-offs)
        ]  # don't include
        trigger_idx = np.array(trigger_idx)
        return trigger_idx, trigger_len
    else:
        return trigger_idx, 0


def calc_trigger_average(neural_data, trigger_idx, trigger_len, pre_event, post_event):
    num_frame = neural_data.shape[0]
    event_trigger_data = []
    for idx in trigger_idx:
        event_start = idx - abs(pre_event)  # use abs to make sure neg value
        event_end = idx + trigger_len + post_event
        if event_end <= num_frame and event_start >= 0:
            event_trigger_data.append(neural_data[event_start:event_end])
    # Convert to numpy array (num_event, cell_number, event_time)
    event_trigger_data = np.array(event_trigger_data)
    return event_trigger_data


def ETA(
    neural_data: FluoData,
    behaviors_data: BehaviorData,
    output_dir: str,
    iscell: IscellData = None,
    params: dict = None,
    **kwargs,
) -> dict(mean=TimeSeriesData):
    function_id = output_dir.split("/")[-1]
    print("start ETA:", function_id)

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

    Y = Y[:, params["event_col_index"]]

    # calculate Triggers
    [trigger_idxs, trigger_len] = calc_trigger(
        Y, params["trigger_type"], params["trigger_threshold"]
    )

    # calculate Triggered average
    event_trigger_data = calc_trigger_average(
        X, trigger_idxs, trigger_len, params["pre_event"], params["post_event"]
    )

    # (cell_number, event_time_lambda)
    if len(event_trigger_data) > 0:
        mean = np.mean(event_trigger_data, axis=0)
        std = np.std(event_trigger_data, axis=0)
        sem = np.std(event_trigger_data, axis=0) / np.sqrt(len(event_trigger_data))
        mean = mean.transpose()
        std = std.transpose()
        sem = sem.transpose()
    else:
        assert False, "Output data size is 0"

    nwbfile = {}
    nwbfile[NWBDATASET.POSTPROCESS] = {
        function_id: {
            "mean": mean,
            "std": std,
            "sem": sem,
            "num_sample": [len(mean)],
        }
    }
    min_value = np.min(mean, axis=1, keepdims=True)
    max_value = np.max(mean, axis=1, keepdims=True)
    norm_mean = (mean - min_value) / (max_value - min_value)

    info = {}
    info["mean"] = TimeSeriesData(
        mean,
        std=std,
        sem=sem,
        index=list(range(params["pre_event"], params["post_event"] + trigger_len)),
        cell_numbers=cell_numbers if iscell is not None else None,
        file_name="mean",
    )
    info["mean_heatmap"] = HeatMapData(
        norm_mean,
        columns=list(
            np.arange(params["pre_event"], params["post_event"] + trigger_len)
        ),
        file_name="mean_heatmap",
    )
    info["nwbfile"] = nwbfile
    return info
