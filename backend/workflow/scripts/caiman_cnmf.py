import sys
sys.path.append('../../../optinist')
import pickle

from wrappers.data_wrapper import *
from wrappers.args_check import args_check
from wrappers.caiman_wrapper import caiman_cnmf



if __name__ == "__main__":

    with open(snakemake.input[0], "rb") as f:
        input_file = pickle.load(f)

    image_data = input_file.values()
    params = snakemake.config["rules"]["caiman_cnmf"]["param"]


    info = caiman_cnmf(*image_data)

    # output to pickle
    with open(snakemake.output[0], 'wb') as f:
        pickle.dump(info, f)
    