from optinist.api.dir_path import DIRPATH
from optinist.api.experiment.experiment import ExptConfig, ExptFunction
from optinist.api.experiment.experiment_writer import ExptConfigWriter
from optinist.api.utils.filepath_creater import create_filepath, join_filepath
from optinist.api.workflow.workflow import Node, NodeData, RunItem

expt_filepath = create_filepath(
    join_filepath([DIRPATH.ROOT_DIR, "test_data"]),
    DIRPATH.EXPERIMENT_YML
)


node_data = NodeData(
    label="a",
    param={},
    path="",
    type=""
)


nodeDict = {
    "test1": Node(
        id="node_id",
        type="a",
        data=node_data,
        position=None,
        style=None,
    )
}


def test_config() -> ExptConfig:
    runItem = RunItem(
        name="New Flow",
        nodeDict={},
        edgeDict={},
        snakemakeParam={},
        nwbParam={},
        forceRunList=[],
    )

    expt_config = ExptConfigWriter.config(
        expt_filepath,
        runItem.name,
        runItem.nodeDict,
        runItem.edgeDict,
    )

    assert isinstance(expt_config, ExptConfig)
    assert isinstance(expt_config.function, dict)
    assert len(expt_config.function) == 0

    return expt_config


def test_add_run_info():
    expt_config = ExptConfigWriter.add_run_info(
        expt_config=test_config(),
        nodeDict=nodeDict,
        edgeDict=None,
    )

    assert len(expt_config.nodeDict) == 1


def test_function_from_nodeDict():
    expt_function = ExptConfigWriter.function_from_nodeDict(nodeDict)

    assert isinstance(expt_function, dict)
    assert isinstance(expt_function["node_id"], ExptFunction)
