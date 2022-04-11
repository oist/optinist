from optinist.api.dir_path import DIRPATH
from optinist.api.experiment.experiment import ExptConfig
from optinist.api.experiment.experiment_writer import ExptConfigWriter
from optinist.api.utils.filepath_creater import create_filepath, join_filepath
from optinist.api.workflow.workflow import RunItem


def test_filepath():
    exp_filepath = create_filepath(
        join_filepath([DIRPATH.ROOT_DIR, "test_data"]),
        DIRPATH.EXPERIMENT_YML
    )
    return exp_filepath


def test_create():
    exp_filepath = test_filepath()

    runItem = RunItem(
        name="New Flow",
        nodeList=[],
        edgeList=[],
        snakemakeParam={},
        nwbParam={},
        forceRunList=[],
    )

    exp_config = ExptConfigWriter.create(
        exp_filepath,
        runItem.name,
        runItem.nodeList,
        runItem.edgeList,
    )

    assert isinstance(exp_config, ExptConfig)
    assert isinstance(exp_config.function, dict)
    assert len(exp_config.function) == 0