from studio.core.snakemake.snakemake_rule import SmkRule
from studio.core.workflow.workflow import Edge, Node, NodeData, NodePosition

unique_id = "test"

node = Node(
    id="input_0",
    type="type",
    data=NodeData(
        label="label",
        param={},
        path="path",
        type="type",
        fileType="fileType",
        hdf5Path="hdf5Path",
    ),
    position=NodePosition(
        x=0,
        y=0,
    ),
    style={},
)

nodeDict = {"node1": node}

edgeDict = {
    "edge1": Edge(
        id="id",
        type="type",
        animated=False,
        source="input_0",
        sourceHandle="input_0--image--ImageData",
        target="suite2p",
        targetHandle="suite2p--image--ImageData",
        style={},
    ),
}


def test_SmkSetfile_image():
    rule = SmkRule(
        unique_id=unique_id,
        node=node,
        edgeDict=edgeDict,
        nwbfile={},
    ).image()

    assert rule.type == "image"


def test_SmkSetfile_csv():
    rule = SmkRule(
        unique_id=unique_id,
        node=node,
        edgeDict=edgeDict,
        nwbfile={},
    ).csv()
    assert rule.type == "csv"


def test_SmkSetfile_hdf5():
    rule = SmkRule(
        unique_id=unique_id,
        node=node,
        edgeDict=edgeDict,
        nwbfile={},
    ).hdf5()

    assert rule.type == "hdf5"


def test_SmkSetfile_algo():
    rule = SmkRule(
        unique_id=unique_id,
        node=node,
        edgeDict=edgeDict,
    ).algo(nodeDict=nodeDict)

    assert rule.type == node.data.label
    assert rule.path == node.data.path
