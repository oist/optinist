from dataclasses import dataclass


@dataclass
class NodeType:
    IMAGE: str = "ImageFileNode"
    CSV: str = "CsvFileNode"
    FLUO: str = "FluoFileNode"
    BEHAVIOR: str = "BehaviorFileNode"
    HDF5: str = "HDF5FileNode"
    ALGO: str = "AlgorithmNode"