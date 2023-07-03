from pydantic import BaseModel


class SnakemakeParams(BaseModel):
    use_conda: bool
    cores: int
    forceall: bool
    forcetargets: bool
    lock: bool
