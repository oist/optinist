from pydantic import BaseModel


class SnakemakeParams(BaseModel):
    use_conda: bool
    cores: int
    forceall: bool
    # These are currently (v2.3.0) not used, but kept for future compatibility
    # forcetargets: bool
    # lock: bool
