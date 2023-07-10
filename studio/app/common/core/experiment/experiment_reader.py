from typing import Dict

import yaml

from studio.app.common.core.experiment.experiment import ExptConfig, ExptFunction
from studio.app.common.core.workflow.workflow import (
    Edge,
    Node,
    NodeData,
    NodePosition,
    OutputPath,
    Style,
)


class ExptConfigReader:
    @classmethod
    def read(cls, filepath) -> ExptConfig:
        with open(filepath, "r") as f:
            config = yaml.safe_load(f)

        return ExptConfig(
            workspace_id=config["workspace_id"],
            unique_id=config["unique_id"],
            name=config["name"],
            started_at=config["started_at"],
            finished_at=config.get("finished_at"),
            success=config.get("success", "running"),
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
                started_at=value.get("started_at"),
                finished_at=value.get("finished_at"),
                success=value.get("success", "running"),
                hasNWB=value["hasNWB"],
                message=value.get("message"),
                outputPaths=cls.read_output_paths(value.get("outputPaths")),
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
    def read_output_paths(cls, config) -> Dict[str, OutputPath]:
        if config:
            return {
                key: OutputPath(
                    path=value["path"],
                    type=value["type"],
                    max_index=value["max_index"],
                )
                for key, value in config.items()
            }
        else:
            return None

    @classmethod
    def rename(cls, filepath, new_name: str) -> ExptConfig:
        with open(filepath, "r") as f:
            config = yaml.safe_load(f)
            config["name"] = new_name

        with open(filepath, "w") as f:
            yaml.dump(config, f)

        return ExptConfig(
            workspace_id=config["workspace_id"],
            unique_id=config["unique_id"],
            name=config["name"],
            started_at=config.get("started_at"),
            finished_at=config.get("finished_at"),
            success=config.get("success", "running"),
            hasNWB=config["hasNWB"],
            function=cls.read_function(config["function"]),
            nodeDict=cls.read_nodeDict(config["nodeDict"]),
            edgeDict=cls.read_edgeDict(config["edgeDict"]),
        )
