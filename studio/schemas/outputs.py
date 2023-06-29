from dataclasses import dataclass
from typing import Dict, List, Union


@dataclass
class OutputData:
    data: Union[List, Dict, str]
    columns: List[str] = None
    index: List[str] = None


@dataclass
class JsonTimeSeriesData(OutputData):
    xrange: list = None
    std: Dict[str, dict] = None
