import traceback
import gc
import copy

from optinist.api.snakemake.smk import Rule
from optinist.api.utils.filepath_creater import create_directory, join_filepath
from optinist.api.pickle.pickle_reader import PickleReader
from optinist.api.pickle.pickle_writer import PickleWriter
from optinist.api.nwb.nwb_creater import merge_nwbfile, save_nwb
from optinist.wrappers import wrapper_dict


class Runner:
    @classmethod
    def run(cls, __rule_config: Rule, last_output):
        try:
            input_info = _read_input_info(__rule_config.input)

            _change_dict_key_exist(input_info, __rule_config)

            for key in list(input_info):
                if key != "nwbfile" and key not in __rule_config.return_arg.values():
                    input_info.pop(key)

            output_info = _execute_function(
                __rule_config.path,
                __rule_config.params,
                input_info
            )

            # ファイル保存先
            output_dir = join_filepath(__rule_config.output.split("/")[:-1])
            create_directory(output_dir)

            # nwbfileの設定
            if "nwbfile" not in output_info:
                output_info["nwbfile"] = input_info["nwbfile"]

            # NWB保存
            if __rule_config.output in last_output:
                output_dir = __rule_config.output.split(".")[0]
                save_nwb(output_info['nwbfile'], output_dir)

            # 結果を保存
            PickleWriter.write(__rule_config.output, output_info)

            print("output: ", __rule_config.output)

            del input_info, output_info
            gc.collect()

        except Exception as e:
            error_message  = list(traceback.TracebackException.from_exception(e).format())[-2:]
            PickleWriter.write(__rule_config.output, error_message)


def _execute_function(path, params, input_info):
    wrapper = _dict2leaf(
        wrapper_dict,
        path.split('/')
    )
    func = copy.deepcopy(wrapper["function"])
    output_info = func(params=params, **input_info)
    del func
    gc.collect()
    return output_info


def _change_dict_key_exist(input_info, rule_config: Rule):
    for return_name, arg_name in rule_config.return_arg.items():
        if return_name in input_info:
            input_info[arg_name] = input_info.pop(return_name)


def _read_input_info(input_files):
    input_info = {}
    for filepath in input_files:
        data = PickleReader.read(filepath)

        input_info = dict(list(data.items()) + list(input_info.items()))
        if 'nwbfile' in input_info:
            input_info['nwbfile'] = merge_nwbfile(
                input_info['nwbfile'],
                data['nwbfile']
            )

    return input_info


def _dict2leaf(root_dict: dict, path_list):
    path = path_list.pop(0)
    if len(path_list) > 0:
        return _dict2leaf(root_dict[path], path_list)
    else:
        return root_dict[path]
