from .dimension_reduction import dimension_reduction_wrapper_dict
from .neural_population_analysis import neural_population_analysis_wrapper_dict
from .neural_decoding import neural_decoding_wrapper_dict


optinist_wrapper_dict = {
    'optinist': {
        'dimension_reduction': dimension_reduction_wrapper_dict ,
        'neural_population_analysis': neural_population_analysis_wrapper_dict,
        'neural_decoding': neural_decoding_wrapper_dict,
    }
}