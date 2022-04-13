from optinist.api.dir_path import DIRPATH
from optinist.api.snakemake.smk import Rule


class SmkConfigReader:
    @classmethod
    def read(cls, rule):
        return Rule(
            rule_file=rule["rule_file"],
            input=rule["input"],
            return_arg=rule["return_arg"],
            params=rule["params"],
            output=rule["output"],
            type=rule["type"],
            nwbfile=rule["nwbfile"],
            hdf5Path=rule["hdf5Path"],
            path=rule["path"],
        )
