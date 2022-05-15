import yaml

from optinist.api.experiment.experiment import (
    ExptConfig,
    ExptFunction,
)
from optinist.api.workflow.workflow import (
    Edge,
    Node,
    NodeData,
    NodePosition,
    Style
)


class ExptConfigReader:
    @classmethod
    def read(cls, filepath) -> ExptConfig:
        with open(filepath, "r") as f:
            config = yaml.safe_load(f)

        return ExptConfig(
            unique_id=config["unique_id"],
            name=config["name"],
            timestamp=config["timestamp"],
            function=cls.function_read(config["function"]),
            nodeDict=cls.nodeDict_read(config["nodeDict"]),
            edgeDict=cls.edgeDict_read(config["edgeDict"]),
        )

    @classmethod
    def function_read(cls, config) -> ExptFunction:
        return {
            key: ExptFunction(
                unique_id=value["unique_id"],
                name=value["name"],
                success=value["success"],
            )
            for key, value in config.items()
        }

    @classmethod
    def nodeDict_read(cls, config) -> Node:
        return { 
            key: 
            Node(
                id=key,
                type=value["type"],
                data=NodeData(**value["data"]),
                position=NodePosition(**value["position"]),
                style=Style(**value["style"])
            )
            for key, value in config.items()
        }

    @classmethod
    def edgeDict_read(cls, config) -> Edge:
        return {
            key:
            Edge(
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
