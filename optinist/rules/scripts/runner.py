import os
import traceback
import gc
import copy

from optinist.wrappers import wrapper_dict
from optinist.api.snakemake.smk import Rule
from optinist.api.pickle.pickle_reader import PickleReader
from optinist.api.pickle.pickle_writer import PickleWriter
from optinist.api.utils.filepath_creater import join_filepath
from optinist.api.nwb.nwb_creater import merge_nwbfile, save_nwb


class Runner:
    @classmethod
    def run(cls, __rule: Rule, last_output):
        try:
            input_info = _read_input_info(__rule.input)

            _change_dict_key_exist(input_info, __rule)

            nwbfile = input_info['nwbfile']

            # input_info
            for key in list(input_info):
                if key not in __rule.return_arg.values():
                    input_info.pop(key)

            # output_info
            output_info = _execute_function(
                __rule.path,
                __rule.params,
                input_info
            )

            # nwbfileの設定
            output_info['nwbfile'] = _save_func_nwb(
                f"{__rule.output.split('.')[0]}.nwb",
                __rule.type,
                nwbfile,
                output_info,
            )

            # 各関数での結果を保存
            PickleWriter.write(__rule.output, output_info)

            # NWB全体保存
            if __rule.output in last_output:
                # 全体の結果を保存する
                path = join_filepath(os.path.dirname(os.path.dirname(__rule.output)))
                path = join_filepath([path, f"all_{__rule.type}.nwb"])
                _save_all_nwb(
                    path,
                    output_info['nwbfile']
                )

            print("output: ", __rule.output)

            del input_info, output_info
            gc.collect()

        except Exception as e:
            PickleWriter.write(
                __rule.output,
                list(traceback.TracebackException.from_exception(e).format())[-2:]
            )


def _save_func_nwb(save_path, name, nwbfile, output_info):
    if "nwbfile" in output_info:
        nwbfile[name] = output_info["nwbfile"]
        save_nwb(
            save_path,
            nwbfile["input"],
            output_info["nwbfile"],
        )
    return nwbfile


def _save_all_nwb(save_path, all_nwbfile):
    input_nwbfile = all_nwbfile["input"]
    all_nwbfile.pop("input")
    nwbfile = {}
    for x in all_nwbfile.values():
        nwbfile = merge_nwbfile(nwbfile, x)
    save_nwb(
        save_path,
        input_nwbfile,
        nwbfile
    )


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
        load_data = PickleReader.read(filepath)
        input_info = dict(list(load_data.items()) + list(input_info.items()))
    return input_info


def _dict2leaf(root_dict: dict, path_list):
    path = path_list.pop(0)
    if len(path_list) > 0:
        return _dict2leaf(root_dict[path], path_list)
    else:
        return root_dict[path]
