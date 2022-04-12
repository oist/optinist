
from dataclasses import dataclass
from typing import Dict, Union

from pydantic import BaseModel


@dataclass
class Rule:
    input: list
    output: str
    type: str
    rule_file: str
    return_arg: str
    params: list
    nwbfile: list


@dataclass
class SnakemakeConfig:
    last_output: list
    rules: Rule


@dataclass
class FlowConfig:
    rules: Dict[str, Rule]
    last_output: list


class ForceRun(BaseModel):
    nodeId: str
    name: str


@dataclass
class Rule:
    rule_file: str
    input: str
    return_arg: Union[str, Dict[str, str]]
    params: dict
    output: str
    type: str
    nwbfile: dict = None
    hdf5Path: str = None
    path: str = None