import pytest
import os
from dataclasses import asdict

from optinist.api.config.config_reader import ConfigReader
from optinist.api.config.config_writer import ConfigWriter
from optinist.api.dir_path import DIRPATH
from optinist.api.snakemake.smk import FlowConfig, SmkParam
from optinist.api.snakemake.snakemake_executor import create_edge_dict, delete_depemdemcies, get_dependencies_graph
from optinist.api.snakemake.snakemake_reader import RuleConfigReader


def test_get_params():
    return SmkParam(
        use_conda=True,
        cores=2,
        forceall=False,
        forcetargets=True,
        lock=False,
    )


def test_get_dependencies_graph():
    # configファイルをtest_data→ROOT_DIRに送る
    test_dirpath = os.path.join(DIRPATH.ROOT_DIR, 'test_data')
    smk_dirpath = os.path.join(test_dirpath, 'snakemake')
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
                    input_list.append(os.path.join(smk_dirpath, x))
            else:
                input_list = os.path.join(smk_dirpath, rule[key_name])

            rule[key_name] = input_list

        new_rules[key] = RuleConfigReader.read(rule)

    new_config = FlowConfig(
        last_output=[
            os.path.join(smk_dirpath, x)
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

    params = test_get_params()

    edges, file_graph = get_dependencies_graph(params)

    # print(edges)
    # for x in edges:
    #     print(x[0], " → ", x[1])
    #     print(file_graph[x[0]], " → ", file_graph[x[1]])

    assert edges
    assert file_graph

    return edges, file_graph


def test_create_edges_dict():
    edges, _ = test_get_dependencies_graph()
    edge_dict = create_edge_dict(edges)

    assert isinstance(edge_dict, dict)
    assert isinstance(edge_dict[1], list)


def test_create_dummy_file(dirname, del_item):
    test_dirpath = os.path.join(DIRPATH.ROOT_DIR, 'test_data')
    smk_dirpath = os.path.join(test_dirpath, 'snakemake')
    dir_filepath = os.path.join(smk_dirpath, dirname)
    os.makedirs(dir_filepath, exist_ok=True)

    del_filepath = os.path.join(dir_filepath, del_item)
    with open(del_filepath, "w") as f:
        f.write("")

    return del_filepath


def test_delete_dependencies():
    "{1: [0], 2: [1], 3: [2]}"
    edges, file_graph = test_get_dependencies_graph()
    edge_dict = create_edge_dict(edges)

    del_filepath = test_create_dummy_file("2", "suite2p_roi.pkl")
    assert os.path.exists(del_filepath)

    delete_depemdemcies(edge_dict, file_graph, del_filepath)
    assert not os.path.exists(del_filepath)

    del_list = []
    del_filepath = test_create_dummy_file("2", "suite2p_roi.pkl")
    del_list.append(del_filepath)
    assert os.path.exists(del_filepath)
    del_filepath = test_create_dummy_file("0", "data_endoscope.pkl")
    del_list.append(del_filepath)
    assert os.path.exists(del_filepath)
    delete_depemdemcies(edge_dict, file_graph, del_filepath)

    for item in del_list:
        assert not os.path.exists(item)
