from optinist.wrappers.optinist.basic_neural_analysis import (
    basic_neural_analysis_wrapper_dict,
)
from optinist.wrappers.optinist.dimension_reduction import (
    dimension_reduction_wrapper_dict,
)
from optinist.wrappers.optinist.neural_decoding import neural_decoding_wrapper_dict
from optinist.wrappers.optinist.neural_population_analysis import (
    neural_population_analysis_wrapper_dict,
)

optinist_wrapper_dict = {
    "optinist": {
        "basic_neural_analysis": basic_neural_analysis_wrapper_dict,
        "dimension_reduction": dimension_reduction_wrapper_dict,
        "neural_population_analysis": neural_population_analysis_wrapper_dict,
        "neural_decoding": neural_decoding_wrapper_dict,
    }
}
