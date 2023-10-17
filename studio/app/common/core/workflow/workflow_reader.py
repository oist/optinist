from typing import Dict

import yaml

from studio.app.common.core.workflow.workflow import (
    Edge,
    Node,
    NodeData,
    NodePosition,
    Style,
)
from studio.app.common.schemas.workflow import WorkflowConfig


class WorkflowConfigReader:
    @classmethod
    def read(cls, filepath) -> WorkflowConfig:
        with open(filepath, "r") as f:
            config = yaml.safe_load(f)

        return WorkflowConfig(
            nodeDict=cls.read_nodeDict(config["nodeDict"]),
            edgeDict=cls.read_edgeDict(config["edgeDict"]),
        )

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
