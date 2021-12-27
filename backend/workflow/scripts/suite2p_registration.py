import sys
sys.path.append('../../../optinist')
import pickle
from wrappers.data_wrapper import *
from wrappers.suite2p_wrapper import suite2p_registration


if __name__ == "__main__":
    with open(snakemake.input[0], "rb") as f:
        input_file = pickle.load(f)

    input_data = input_file.values()
    print(input_data)
    params = snakemake.config["rules"]["suite2p_registration"]["param"]
    
    info = suite2p_registration(*input_data, params=params)

    # output to pickle
    with open(snakemake.output[0], 'wb') as f:
        pickle.dump(info, f)