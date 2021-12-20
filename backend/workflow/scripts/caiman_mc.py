import sys
sys.path.append('../../../optinist')
import pickle
from wrappers.data_wrapper import *
from wrappers.caiman_wrapper import caiman_mc



if __name__ == "__main__":
    input_file = snakemake.input[0]
    input_imagedata = ImageData(input_file, '')
    output_file = snakemake.output[0]
    params = snakemake.config["rules"]["caiman_mc"]["param"]
    info = caiman_mc(input_imagedata, params=params)
    
    # output to pickle
    with open(output_file, 'wb') as f:
        pickle.dump(info, f)

