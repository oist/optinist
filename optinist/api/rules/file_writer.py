import h5py

from optinist.api.dataclass.dataclass import CsvData, ImageData, TimeSeriesData
from optinist.api.nwb.nwb import NWBDATASET
from optinist.api.snakemake.smk import Rule
from optinist.routers.model import FILETYPE


class FileWriter:
    @classmethod
    def csv(cls, rule_config: Rule, nodeType):
        info = {
            rule_config.return_arg: CsvData(rule_config.input, rule_config.params, "")
        }
        nwbfile = rule_config.nwbfile

        if nodeType == FILETYPE.CSV:
            if NWBDATASET.TIMESERIES not in nwbfile:
                nwbfile[NWBDATASET.TIMESERIES] = {}
            nwbfile[NWBDATASET.TIMESERIES][rule_config.return_arg] = info[
                rule_config.return_arg
            ]
        elif nodeType == FILETYPE.BEHAVIOR:
            if NWBDATASET.BEHAVIOR not in nwbfile:
                nwbfile[NWBDATASET.BEHAVIOR] = {}
            nwbfile[NWBDATASET.BEHAVIOR][rule_config.return_arg] = info[
                rule_config.return_arg
            ]
        else:
            assert False, "NodeType doesn't exsits"

        nwbfile.pop("image_series", None)
        info["nwbfile"] = {"input": nwbfile}
        return info

    @classmethod
    def image(cls, rule_config: Rule):
        info = {rule_config.return_arg: ImageData(rule_config.input, "")}
        nwbfile = rule_config.nwbfile
        nwbfile["image_series"]["external_file"] = info[rule_config.return_arg]
        info["nwbfile"] = {"input": nwbfile}
        return info

    @classmethod
    def hdf5(cls, rule_config: Rule):
        nwbfile = rule_config.nwbfile

        with h5py.File(rule_config.input, "r") as f:
            data = f[rule_config.hdf5Path][:]

        if data.ndim == 3:
            info = {rule_config.return_arg: ImageData(data, "")}
            nwbfile["image_series"]["external_file"] = info[rule_config.return_arg]
            info["nwbfile"] = {}
            info["nwbfile"][FILETYPE.IMAGE] = nwbfile
        elif data.ndim == 2:
            info = {rule_config.return_arg: TimeSeriesData(data, "")}

            if NWBDATASET.TIMESERIES not in nwbfile:
                nwbfile[NWBDATASET.TIMESERIES] = {}

            nwbfile[NWBDATASET.TIMESERIES][rule_config.return_arg] = info[
                rule_config.return_arg
            ]
            nwbfile.pop("image_series", None)
            info["nwbfile"] = {"input": nwbfile}

        return info
