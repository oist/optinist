
from dataclasses import asdict, dataclass, field
import os
from typing import Dict, List, Union

from pydantic import BaseModel

from optinist.api.dir_path import DIRPATH


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

    def to_cluster_rule(self):
        # 非クラスタマシン用のRuleインスタンスをもとに、
        # 新たにクラスタマシン用のRuleを作成する

        # inputとoutputのpathをクラスタマシン用に変換する
        cluster_input = [self.replace_upper_directory(input) for input in self.input]
        cluster_output = self.replace_upper_directory(self.output)

        cluster_rule_arg_dict = {
            **asdict(self),
            "input": cluster_input,
            "output": cluster_output,
        }
        cluster_rule = Rule(
            **cluster_rule_arg_dict,
        )
        return cluster_rule

    @classmethod
    def replace_upper_directory(_cls, path):
        # 非クラスタマシン用に設定されたディレクトリを
        # クラスタ用に設定されたディレクトリに書き換える
        acc = []
        input_dir_checker = os.path.basename(DIRPATH.INPUT_DIR)
        output_dir_checker = os.path.basename(DIRPATH.OUTPUT_DIR)
        remaining = path

        while True:
            remaining, tmp = os.path.split(remaining)
            to_be_checked = os.path.basename(remaining)
            acc = [tmp] + acc
            if to_be_checked == input_dir_checker:
                result = os.path.join(DIRPATH.CLUSTER_INPUT_DIR, *acc)
                break
            elif to_be_checked == output_dir_checker:
                result = os.path.join(DIRPATH.CLUSTER_OUTPUT_DIR, *acc)
                break

        return result


@dataclass
class FlowConfig:
    rules: Dict[str, Rule]
    last_output: list


class ForceRun(BaseModel):
    nodeId: str
    name: str


@dataclass
class SmkParam:
    use_conda: bool
    cores: int
    forceall: bool
    forcetargets: bool
    lock: bool
    forcerun: List[str] = field(default_factory=list)
