import os
import yaml
from datetime import datetime
from dataclasses import dataclass, asdict

from optinist.cui_api.dir_path import DIRPATH
from optinist.cui_api.filepath_creater import join_filepath
from optinist.cui_api.config_writer import ConfigWriter

from typing import Dict, List

@dataclass
class ExpFunction:
    unique_id: str
    name: str
    success: str

@dataclass
class NodeData:
    label: str
    param: dict
    path: list or str
    type: str
    fileType: str = None

@dataclass
class NodePosition:
    x: int
    y: int

@dataclass
class Style:
    border: str = None
    height: int = None
    padding: int = None
    width: int = None
    borderRadius: int = None

@dataclass
class Node:
    id: str
    type: str
    data: NodeData
    position: NodePosition
    style: Style

@dataclass
class Edge:
    id: str
    type: str
    animated: bool
    source: str
    sourceHandle: str
    target: str
    targetHandle: str
    style: Style

@dataclass
class ExpConfig:
    timestamp: str
    name: str
    unique_id: str
    function: Dict[str, ExpFunction]
    nodeList: List[Node]
    edgeList: List[Edge]


class ExpConfigReader:
    @classmethod
    def read(cls, filepath) -> ExpConfig:
        with open(filepath, "r") as f:
            config = yaml.safe_load(f)

        return ExpConfig(
            timestamp=config["timestamp"],
            name=config["name"],
            unique_id=config["unique_id"],
            function=cls.function_read(config["function"]),
            nodeList=cls.nodeList_read(config["nodeList"]),
            edgeList=cls.edgeList_read(config["edgeList"]),
        )

    @classmethod
    def function_read(cls, config) -> ExpFunction:
        return {
            key: ExpFunction(
                unique_id=value["unique_id"],
                name=value["name"],
                success=value["success"],
            )
            for key, value in config.items()
        }

    @classmethod
    def nodeList_read(cls, config) -> Node:
        return [
            Node(
                id=value["id"],
                type=value["type"],
                data=NodeData(**value["data"]),
                position=NodePosition(**value["position"]),
                style=Style(**value["style"])
            )
            for value in config
        ]
        

    @classmethod
    def edgeList_read(cls, config) -> Edge:
        return [
            Edge(
                id=value["id"],
                type=value["type"],
                animated=value["animated"],
                source=value["source"],
                sourceHandle=value["sourceHandle"],
                target=value["target"],
                targetHandle=value["targetHandle"],
                style=Style(**value["style"]),
            )
            for value in config
        ]


def create_function_from_nodeList(nodeList):
    return {
        node["id"]: ExpFunction(
            unique_id=node["id"],
            name=node["data"]["label"],
            success="success" if node["data"]["type"] == "input" else "running",
        )
        for node in nodeList
    }


def create_exp_config(unique_id, runItem):
    return ExpConfig(
        unique_id=unique_id,
        name=runItem.name,
        timestamp=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        nodeList=ExpConfigReader.nodeList_read(runItem.nodeList),
        edgeList=ExpConfigReader.edgeList_read(runItem.edgeList),
        function=create_function_from_nodeList(runItem.nodeList),
    )


def add_run_info(exp_config, runItem):
    exp_config.timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # 関数を追加の可能性
    exp_config.nodeList += ExpConfigReader.nodeList_read(runItem.nodeList)
    exp_config.edgeList += ExpConfigReader.edgeList_read(runItem.edgeList)

    for node in runItem.nodeList:
        exp_config.function[node["id"]] = ExpFunction(
            unique_id=node["id"],
            name=node["data"]["label"],
            success="success" if node["data"]["type"] == "input" else "running",
        )
    return exp_config


def save_exp_config(filepath, exp_config):
    with open(filepath, "w") as f:
        yaml.dump(exp_config, f)


def exp_config_writer(unique_id, rule_config, runItem):
    exp_filepath = join_filepath([DIRPATH.BASE_DIR, unique_id, "experiment.yaml"])
    if os.path.exists(exp_filepath):
        exp_config = ExpConfigReader.read(exp_filepath)
        exp_config = add_run_info(exp_config, runItem)
    else:
        exp_config = create_exp_config(unique_id, runItem)

    ConfigWriter.write(
        dirname=join_filepath([DIRPATH.BASE_DIR, unique_id]),
        filename="experiment.yaml",
        config=asdict(exp_config),
    )
