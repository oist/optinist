from optinist.api.dir_path import DIRPATH
from optinist.api.experiment.experiment import ExptConfig
from optinist.api.experiment.experiment_writer import _create_config
from optinist.api.utils.filepath_creater import create_filepath, join_filepath
from optinist.api.workflow.workflow import RunItem


def test_filepath():
    exp_filepath = create_filepath(
        join_filepath([DIRPATH.ROOT_DIR, "test_data"]),
        DIRPATH.EXPERIMENT_YML
    )
    return exp_filepath


def test_create_config():
    exp_filepath = test_filepath()

    runItem = RunItem(
        name="New Flow",
        nodeDict={},
        edgeDict={},
        snakemakeParam={},
        nwbParam={},
        forceRunList=[],
    )

    exp_config = _create_config(
        exp_filepath,
        runItem.name,
        runItem.nodeDict,
        runItem.edgeDict,
    )

    assert isinstance(exp_config, ExptConfig)
    assert isinstance(exp_config.function, dict)
    assert len(exp_config.function) == 0