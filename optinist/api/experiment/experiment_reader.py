from typing import Dict

import yaml

from optinist.api.experiment.experiment import ExptConfig, ExptFunction
from optinist.api.workflow.workflow import Edge, Node, NodeData, NodePosition, Style


class ExptConfigReader:
    @classmethod
    def read(cls, filepath) -> ExptConfig:
        with open(filepath, "r") as f:
            config = yaml.safe_load(f)

        return ExptConfig(
            unique_id=config["unique_id"],
            name=config["name"],
            timestamp=config["timestamp"],
            hasNWB=config["hasNWB"],
            function=cls.read_function(config["function"]),
            nodeDict=cls.read_nodeDict(config["nodeDict"]),
            edgeDict=cls.read_edgeDict(config["edgeDict"]),
        )

    @classmethod
    def read_function(cls, config) -> Dict[str, ExptFunction]:
        return {
            key: ExptFunction(
                unique_id=value["unique_id"],
                name=value["name"],
                success=value["success"],
                hasNWB=value["hasNWB"],
            )
            for key, value in config.items()
        }

    @classmethod
    def read_nodeDict(cls, config) -> Dict[str, Node]:
        return {
            key: Node(
                id=key,
                type=value["type"],
                data=NodeData(**value["data"]),
                position=NodePosition(**value["position"]),
                style=Style(**value["style"]),
            )
            for key, value in config.items()
        }

    @classmethod
    def read_edgeDict(cls, config) -> Dict[str, Edge]:
        return {
            key: Edge(
                id=key,
                type=value["type"],
                animated=value["animated"],
                source=value["source"],
                sourceHandle=value["sourceHandle"],
                target=value["target"],
                targetHandle=value["targetHandle"],
                style=Style(**value["style"]),
            )
            for key, value in config.items()
        }

    @classmethod
    def rename(cls, filepath, new_name: str) -> ExptConfig:
        with open(filepath, "r") as f:
            config = yaml.safe_load(f)
            config["name"] = new_name

        with open(filepath, "w") as f:
            yaml.dump(config, f)

        return ExptConfig(
            unique_id=config["unique_id"],
            name=config["name"],
            timestamp=config["timestamp"],
            hasNWB=config["hasNWB"],
            function=cls.read_function(config["function"]),
            nodeDict=cls.read_nodeDict(config["nodeDict"]),
            edgeDict=cls.read_edgeDict(config["edgeDict"]),
        )
