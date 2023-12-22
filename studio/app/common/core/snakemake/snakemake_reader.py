from studio.app.common.core.snakemake.smk import Rule, SmkParam


class RuleConfigReader:
    @classmethod
    def read(cls, rule):
        return Rule(
            input=rule["input"],
            return_arg=rule["return_arg"],
            params=rule["params"],
            output=rule["output"],
            type=rule["type"],
            nwbfile=rule["nwbfile"],
            hdf5Path=rule["hdf5Path"],
            matPath=rule["matPath"],
            path=rule["path"],
        )


class SmkParamReader:
    @classmethod
    def read(cls, params):
        return SmkParam(
            use_conda=params["use_conda"],
            cores=params["cores"],
            forceall=params["forceall"],
            forcetargets=params["forcetargets"],
            lock=params["lock"],
            forcerun=params["forcerun"] if "forcerun" in params else [],
        )
