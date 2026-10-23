from studio.app.common.core.snakemake.smk import Rule, SmkParam
from studio.app.common.core.snakemake.snakemake_reader import (
    RuleConfigReader,
    SmkParamReader,
)

rule_config = {
    "hdf5Path": None,
    "matPath": None,
    "input": ["snakemake/0/data_endoscope.pkl"],
    "nwbfile": None,
    "output": "snakemake/1/suite2p_roi.pkl",
    "params": {
        "batch_size": 500,
        "do_registration": 1,
        "force_sktiff": False,
        "nchannels": 1,
        "nplanes": 1,
    },
    "path": "suite2p/suite2p_file_convert",
    "return_arg": {
        "input_0": "image",
    },
    "type": "suite2p_file_convert",
}

smk_config = {
    "use_conda": False,
    "cores": 2,
    "forceall": True,
    # "forcetargets": True,
    # "lock": False,
}


def test_RuleConfigReader_read():
    output = RuleConfigReader.read(rule_config)

    assert isinstance(output, Rule)


def test_SmkParamReader_read():
    output = SmkParamReader.read(smk_config)

    assert isinstance(output, SmkParam)
