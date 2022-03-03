from wrappers.data_wrapper import *
from wrappers.args_check import args_check
from wrappers.optinist_wrapper.utils import standard_norm


def calc_trigger(behavior_data, trigger_type, trigger_threshold):
    flg = np.array(behavior_data > trigger_threshold, dtype=int)
    if trigger_type == 'up':
        trigger_idx = np.where(np.ediff1d(flg) == 1)[0]
    elif trigger_type == 'down':
        trigger_idx = np.where(np.ediff1d(flg) == -1)[0]
    elif trigger_type == 'cross':
        trigger_idx = np.where(np.ediff1d(flg) != 0)[0]
    else:
        trigger_idx = np.where(np.ediff1d(flg) == 0)[0]

    return trigger_idx


def calc_trigger_average(neural_data, trigger_idx, start_time, end_time):

    num_frame = neural_data.shape[1]

    ind = np.array(range(start_time, end_time+1), dtype = int)

    event_trigger_data = []
    for i in range(len(trigger_idx)):
        target_idx = ind + trigger_idx[i]

        if np.min(target_idx) >= 0 and np.max(target_idx) < num_frame:
            event_trigger_data.append(neural_data[:, target_idx])

    event_trigger_data = np.array(event_trigger_data)

    return event_trigger_data


def ETA(
        neural_data: TimeSeriesData,
        behavior_data: TimeSeriesData,
        iscell: IscellData=None,
        nwbfile: NWBFile=None,
        params: dict=None
    ) -> {}:

    neural_data = neural_data.data
    behavior_data = behavior_data.data

    # data shold be time x component matrix
    if params['transpose_x']:
        neural_data = neural_data.transpose()

    if params['transpose_y']:
        behavior_data = behavior_data.transpose()

    if iscell is not None:
        iscell = iscell.data
        ind  = np.where(iscell > 0)[0]
        neural_data = neural_data[ind, :]
        behavior_data = behavior_data[ind, :]

    behavior_data = behavior_data[:, params['target_index']]

    # calculate Triggers ##################
    trigger_idx = calc_trigger(behavior_data, params['trigger_type'], params['trigger_threshold'])

    # calculate Triggered average ##################
    event_trigger_data = calc_trigger_average(neural_data, trigger_idx, params['start_time'], params['end_time'])
    
    info = {}
    info['event_trigger_data'] = TimeSeriesData(
        event_trigger_data,
        func_name='eta',
        file_name='event_trigger_data'
    )

    return info
