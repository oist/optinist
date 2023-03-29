from dataclasses import dataclass

@dataclass
class ACTION:
    ADD: str = 'add'
    MERGE: str = 'merge'
    DELETE: str = 'delete'