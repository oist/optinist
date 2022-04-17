import pytest
import os
from dataclasses import asdict

from optinist.api.config.config_reader import ConfigReader
from optinist.api.config.config_writer import ConfigWriter
from optinist.api.dir_path import DIRPATH
from optinist.api.snakemake.smk import FlowConfig, SmkParam
from optinist.api.snakemake.snakemake_executor import SmkExecutor, get_dependencies_graph
from optinist.api.snakemake.snakemake_reader import RuleConfigReader

def test_get_dependencies_graph():
    # configファイルをtest_data→ROOT_DIRに送る
    test_dirpath = os.path.join(DIRPATH.ROOT_DIR, 'test_data')
    config_filepath = os.path.join(DIRPATH.ROOT_DIR, 'config.yaml')
    test_config_filepath = os.path.join(test_dirpath, 'config.yaml')

    # configを削除
    if os.path.exists(config_filepath):
        os.remove(config_filepath)
    assert not os.path.exists(config_filepath)

    # test configをROOT_DIRにコピー
    config = ConfigReader.read(test_config_filepath)
    new_rules = {}
    for key, rule in config["rules"].items():
        input_list = []
        for key_name in ["input", "output"]:
            if isinstance(rule[key_name], list):
                for x in rule[key_name]:
                    input_list.append(os.path.join(test_dirpath, "snakemake", x))
            else:
                input_list = os.path.join(test_dirpath, "snakemake", rule[key_name])

            rule[key_name] = input_list

        new_rules[key] = RuleConfigReader.read(rule)

    new_config = FlowConfig(
        last_output=[
            os.path.join(test_dirpath, "snakemake", x)
            for x in config["last_output"]
        ],
        rules=new_rules
    )

    ConfigWriter.write(
        dirname=DIRPATH.ROOT_DIR,
        filename="config.yaml",
        config=asdict(new_config),
    )

    assert os.path.exists(test_config_filepath)

    params = SmkParam(
        use_conda=True,
        cores=2,
        forceall=False,
        forcetargets=True,
        lock=False,
    )
    edges, file_graph = get_dependencies_graph(params)

    assert edges
    assert file_graph
    assert False

    return edges, file_graph


def test_delete_dependencies():
    edges, file_graph = test_get_dependencies_graph()
    assert edges
    get_dependencies_graph(params)
