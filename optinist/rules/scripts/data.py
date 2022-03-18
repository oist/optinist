import sys
sys.path.append('../optinist')

import pickle
import h5py
from wrappers.data_wrapper import ImageData, CsvData, TimeSeriesData
from wrappers.nwb_wrapper.const import NWBDATASET


def save_csv(rule):
    info = {rule["return_arg"]: CsvData(rule["input"], rule["params"], '')}
    nwbfile = rule["nwbfile"]

    if NWBDATASET.TIMESERIES not in nwbfile:
        nwbfile[NWBDATASET.TIMESERIES] = {}

    nwbfile[NWBDATASET.TIMESERIES][rule["return_arg"]] = info[rule["return_arg"]]
    info['nwbfile'] = nwbfile

    save_pickle(rule["output"], info)


def save_image(rule):
    info = {rule["return_arg"]: ImageData(rule["input"], '')}
    nwbfile = rule["nwbfile"]

    # NWB file
    nwbfile['image_series']['external_file'] = info[rule["return_arg"]]
    info['nwbfile'] = nwbfile

    save_pickle(rule["output"], info)


def save_hdf5(rule):
    nwbfile = rule["nwbfile"]

    with h5py.File(rule["input"], "r") as f:
        data = f[rule["hdf5Path"]][:]

    if data.ndim == 3:
        info = {rule["return_arg"]: ImageData(data, '')}
        nwbfile['image_series']['external_file'] = info[rule["return_arg"]]
    elif data.ndim == 2:
        info = {rule["return_arg"]: TimeSeriesData(data, '')}

        if NWBDATASET.TIMESERIES not in nwbfile.keys():
            nwbfile[NWBDATASET.TIMESERIES] = {}
        nwbfile[NWBDATASET.TIMESERIES][rule["return_arg"]] = info[rule["return_arg"]]

    if NWBDATASET.TIMESERIES not in nwbfile:
        nwbfile[NWBDATASET.TIMESERIES] = {}

    nwbfile[NWBDATASET.TIMESERIES][rule["return_arg"]] = info[rule["return_arg"]]
    info['nwbfile'] = nwbfile

    save_pickle(rule["output"], info)


def save_pickle(pickle_path, info):
    with open(pickle_path, 'wb') as f:
        pickle.dump(info, f)


if __name__ == '__main__':
    from rules.utils import run_script
    last_output = snakemake.config["last_output"]

    for rule in snakemake.config["rules"].values():
        if rule["type"] == "csv":
            save_csv(rule)
        elif rule["type"] == "image":
            save_image(rule)
        elif rule["type"] == "hdf5":
            save_hdf5(rule)
