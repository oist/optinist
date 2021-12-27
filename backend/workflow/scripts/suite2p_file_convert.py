import sys
sys.path.append('../../../optinist')
import pickle
from wrappers.data_wrapper import *
from wrappers.suite2p_wrapper import suite2p_file_convert


if __name__ == "__main__":
    input_file = snakemake.input[0]
    input_imagedata = ImageData(input_file, '')
    output_file = snakemake.output[0]
    params = snakemake.config["rules"]["suite2p_file_convert"]["param"]
    info = suite2p_file_convert(input_imagedata, params=params)

    # output to pickle
    with open(output_file, 'wb') as f:
        pickle.dump(info, f)