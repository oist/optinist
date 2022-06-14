# import pytest
# import os
# import shutil

# from optinist.api.dir_path import DIRPATH
# from optinist.api.config.config_reader import ConfigReader
# from optinist.api.config.config_writer import ConfigWriter
# from optinist.api.snakemake.smk import SmkParam
# from optinist.api.snakemake.snakemake_executor import (
#     SmkExecutor,
#     delete_dependencies,
# )
# from optinist.api.utils.filepath_creater import create_directory, join_filepath


# tiff_filename = "test.tif"


# def create_dummy_file(dirname, del_item):
#     relative_filepath = os.path.join('test_data', 'snakemake', dirname, del_item)
#     absolute_dirpath = os.path.join(DIRPATH.OUTPUT_DIR, 'test_data', 'snakemake', dirname)
#     create_directory(absolute_dirpath, delete_dir=True)

#     absolute_filepath = os.path.join(absolute_dirpath, del_item)
#     with open(absolute_filepath, "w") as f:
#         f.write("")

#     return absolute_filepath, relative_filepath


# def create_dummy_params():
#     return SmkParam(
#         use_conda=True,
#         cores=2,
#         forceall=False,
#         forcetargets=True,
#         lock=False,
#     )


# def test_get_dependencies_graph():
#     # configファイルをtest_data→ROOT_DIRに送る
#     test_dirpath = os.path.join(DIRPATH.ROOT_DIR, 'test_data')
#     smk_dirpath = os.path.join(test_dirpath, 'snakemake')

#     config_filepath = os.path.join(DIRPATH.ROOT_DIR, 'config.yaml')
#     test_config_filepath = os.path.join(test_dirpath, 'config.yaml')

#     create_directory(DIRPATH.INPUT_DIR)

#     shutil.copyfile(
#         join_filepath([smk_dirpath, tiff_filename]),
#         join_filepath([DIRPATH.INPUT_DIR, tiff_filename]),
#     )

#     # configを削除
#     if os.path.exists(config_filepath):
#         os.remove(config_filepath)
#     assert not os.path.exists(config_filepath)

#     # # test configをROOT_DIRにコピー
#     config = ConfigReader.read(test_config_filepath)

#     ConfigWriter.write(
#         dirname=DIRPATH.ROOT_DIR,
#         filename="config.yaml",
#         config=config,
#     )

#     assert os.path.exists(test_config_filepath)

#     params = create_dummy_params()

#     smk_executor = SmkExecutor(
#         DIRPATH.SNAKEMAKE_FILEPATH,
#         forceall=params.forceall,
#         cores=params.cores,
#     )
#     edge_dict, file_graph = smk_executor.init_graph()

#     assert edge_dict
#     assert file_graph

#     os.remove(join_filepath([DIRPATH.INPUT_DIR, tiff_filename]))

#     return edge_dict, file_graph


# def test_create_edges_dict():
#     edge_dict, _ = test_get_dependencies_graph()

#     assert isinstance(edge_dict, dict)
#     assert isinstance(edge_dict[1], list)


# def test_delete_dependencies():
#     "{1: [0], 2: [1], 3: [2]}"

#     absolute_filepath, relative_filepath = create_dummy_file("2", "suite2p_roi.pkl")
#     assert os.path.exists(absolute_filepath)

#     params = create_dummy_params()
#     params.forcerun = [relative_filepath]

#     delete_dependencies(params)
#     assert not os.path.exists(absolute_filepath)

#     del_dict = {}
#     absolute_filepath, relative_filepath = create_dummy_file("2", "suite2p_roi.pkl")
#     del_dict[relative_filepath] = absolute_filepath
#     assert os.path.exists(absolute_filepath)

#     absolute_filepath, relative_filepath = create_dummy_file("0", "data_endoscope.pkl")
#     del_dict[relative_filepath] = absolute_filepath
#     assert os.path.exists(absolute_filepath)

#     params = create_dummy_params()
#     params.forcerun = list(del_dict.keys())
#     delete_dependencies(params)

#     for del_item in del_dict.values():
#         assert not os.path.exists(del_item)
