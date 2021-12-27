import sys
sys.path.append('../../../optinist')
import pickle
from wrappers.data_wrapper import *
from wrappers.suite2p_wrapper import suite2p_roi


if __name__ == "__main__":
    with open(snakemake.input[0], "rb") as f:
        input_file = pickle.load(f)

    input_data = input_file.values()
    
    params = snakemake.config["rules"]["suite2p_roi"]["param"]
    info = suite2p_roi(*input_data, params=params)

    output_file = snakemake.output[0]

    # output to pickle
    with open(output_file, 'wb') as f:
        pickle.dump(info, f)