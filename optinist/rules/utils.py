
import traceback
import os
import pickle
from cui_api.utils import join_file_path
import gc
import copy

from wrappers import wrapper_dict
from wrappers.nwb_wrapper import save_nwb, NWBDATASET


def dict2leaf(root_dict: dict, path_list):
    path = path_list.pop(0)
    if len(path_list) > 0:
        return dict2leaf(root_dict[path], path_list)
    else:
        return root_dict[path]


def change_dict_key_exist(d, old_key, new_key):
    if old_key in d:
        d[new_key] = d.pop(old_key)


def merge_nwbfile(old_nwbfile, new_nwbfile):
    for pattern in [
        NWBDATASET.POSTPROCESS,
        NWBDATASET.TIMESERIES,
        NWBDATASET.MOTION_CORRECTION,
        NWBDATASET.ROI,
        NWBDATASET.COLUMN,
        NWBDATASET.FLUORESCENCE,
        NWBDATASET.BEHAVIOR,
    ]:
        if pattern in new_nwbfile:
            if pattern in old_nwbfile:
                old_nwbfile[pattern].update(new_nwbfile[pattern])
            else:
                old_nwbfile[pattern] = new_nwbfile[pattern]
    return old_nwbfile


def get_input_info(input_files):
    input_info = {}
    for path in input_files:
        with open(path, 'rb') as f:
            data = pickle.load(f)
            input_info = dict(list(data.items()) + list(input_info.items()))
            if 'nwbfile' in input_info:
                input_info['nwbfile'] = merge_nwbfile(input_info['nwbfile'], data['nwbfile'])
    return input_info


def run_function(wrapper, params, input_info):
    func = copy.deepcopy(wrapper["function"])
    output_info = func(params=params, **input_info)
    del func
    gc.collect()
    return output_info


def run_script(__func_config, last_output):
    try:
        input_files = __func_config["input"]
        return_arg = __func_config["return_arg"]
        node_type = __func_config["type"]
        params = __func_config["params"]
        wrapper = dict2leaf(wrapper_dict, __func_config["path"].split('/'))
        print(wrapper)

        input_info = get_input_info(input_files)

        for return_name, arg_name in return_arg.items():
            change_dict_key_exist(input_info, return_name, arg_name)

        for key in list(input_info):
            if key != "nwbfile" and key not in return_arg.values():
                input_info.pop(key)

        output_info = run_function(wrapper, params, input_info)

        # ファイル保存先
        output_dir = join_file_path(__func_config["output"].split("/")[:-1])
        os.makedirs(output_dir, exist_ok=True)

        # nwbfileの設定
        if "nwbfile" not in output_info.keys():
            output_info["nwbfile"] = input_info["nwbfile"]

        # NWB保存
        if __func_config["output"] in last_output:
            output_dir = __func_config["output"].split(".")[0]
            save_nwb(output_info['nwbfile'], output_dir)

        # 結果を保存
        with open(__func_config["output"], 'wb') as f:
            pickle.dump(output_info, f)

        print("output: ", __func_config["output"])

        del input_info, output_info
        gc.collect()

    except Exception as e:
        error_message  = list(traceback.TracebackException.from_exception(e).format())[-2:]
        with open(__func_config["output"], 'wb') as f:
            pickle.dump(error_message, f)
