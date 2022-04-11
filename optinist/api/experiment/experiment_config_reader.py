import yaml

from optinist.api.experiment.experiment import (
    ExpConfig,
    ExpFunction,
)
from optinist.api.workflow.workflow import Edge, Node, NodeData, NodePosition, Style


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

