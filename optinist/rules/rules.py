import sys
sys.path.append('../../optinist')
from wrappers.data_wrapper import *
from wrappers.dummy_wrapper.dummy import *
from wrappers import wrapper_dict

rule:
    input:
        "/tmp/optinist/Sue_2x_3000_40_-46/Sue_2x_3000_40_-46.tif"
    output: 
        touch("./files/dummy_image2image/dummy_image2image_out.pkl")
    run:
        info = wrapper_dict["dummy"][config["rules"]["dummy_image2image"]["name"]]["function"](
            ImageData("/tmp/optinist/Sue_2x_3000_40_-46/Sue_2x_3000_40_-46.tif"))
        print(info)

rule:
    input:
        "/tmp/optinist/Sue_2x_3000_40_-46/Sue_2x_3000_40_-46.tif"
    output: 
        touch("./files/dummy_image2time/dummy_image2time_out.pkl")
    run:
        info = wrapper_dict["dummy"][config["rules"]["dummy_image2time"]["name"]]["function"](
            ImageData("/tmp/optinist/Sue_2x_3000_40_-46/Sue_2x_3000_40_-46.tif"))
        print(info)

rule:
    input:
        "./files/dummy_image2image/dummy_image2image_out.pkl"
    output: 
        touch("./files/dummy_image2heat/dummy_image2heat_out.pkl")
    run:
        info = wrapper_dict["dummy"][config["rules"]["dummy_image2heat"]["name"]]["function"](
            ImageData("/tmp/optinist/Sue_2x_3000_40_-46/Sue_2x_3000_40_-46.tif"))
        print(info)