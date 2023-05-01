from optinist.wrappers.optinist_wrapper.basic_neural_analysis import (
    basic_neural_analysis_wrapper_dict,
)
from optinist.wrappers.optinist_wrapper.dimension_reduction import (
    dimension_reduction_wrapper_dict,
)
from optinist.wrappers.optinist_wrapper.neural_decoding import (
    neural_decoding_wrapper_dict,
)
from optinist.wrappers.optinist_wrapper.neural_population_analysis import (
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
