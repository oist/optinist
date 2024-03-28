import h5py

from studio.app.common.core.snakemake.smk import Rule
from studio.app.common.dataclass import CsvData, ImageData, TimeSeriesData
from studio.app.const import FILETYPE
from studio.app.optinist.core.nwb.nwb import NWBDATASET
from studio.app.optinist.dataclass.iscell import IscellData
from studio.app.optinist.dataclass.microscope import MicroscopeData
from studio.app.optinist.routers.mat import MatGetter


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
            assert False, "NodeType doesn't exist"

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

        return cls.get_info_from_array_data(rule_config, nwbfile, data)

    @classmethod
    def mat(cls, rule_config: Rule):
        nwbfile = rule_config.nwbfile
        data = MatGetter.data(rule_config.input, rule_config.matPath)
        return cls.get_info_from_array_data(rule_config, nwbfile, data)

    @classmethod
    def microscope(cls, rule_config: Rule):
        info = {rule_config.return_arg: MicroscopeData(rule_config.input)}
        nwbfile = rule_config.nwbfile
        nwbfile["image_series"]["external_file"] = info[rule_config.return_arg]
        info["nwbfile"] = {"input": nwbfile}
        return info

    @classmethod
    def get_info_from_array_data(cls, rule_config: Rule, nwbfile, data):
        if data.ndim == 3:
            info = {rule_config.return_arg: ImageData(data)}
            nwbfile["image_series"]["external_file"] = info[rule_config.return_arg]
            info["nwbfile"] = {"input": nwbfile}
            info["nwbfile"][FILETYPE.IMAGE] = nwbfile
        elif data.ndim == 2:
            info = {rule_config.return_arg: TimeSeriesData(data)}

            if NWBDATASET.TIMESERIES not in nwbfile:
                nwbfile[NWBDATASET.TIMESERIES] = {}

            nwbfile[NWBDATASET.TIMESERIES][rule_config.return_arg] = info[
                rule_config.return_arg
            ]
            nwbfile.pop("image_series", None)
            info["nwbfile"] = {"input": nwbfile}
        elif data.ndim == 1:
            info = {rule_config.return_arg: IscellData(data)}

            if NWBDATASET.COLUMN not in nwbfile:
                nwbfile[NWBDATASET.COLUMN] = {}

            nwbfile[NWBDATASET.COLUMN][rule_config.return_arg] = info[
                rule_config.return_arg
            ]
            info["nwbfile"] = {"input": nwbfile}
        return info
