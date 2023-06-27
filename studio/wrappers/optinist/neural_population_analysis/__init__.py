from studio.wrappers.optinist.neural_population_analysis.correlation import correlation
from studio.wrappers.optinist.neural_population_analysis.cross_correlation import (
    cross_correlation,
)
from studio.wrappers.optinist.neural_population_analysis.granger import Granger

neural_population_analysis_wrapper_dict = {
    "correlation": {
        "function": correlation,
    },
    "cross_correlation": {
        "function": cross_correlation,
        "conda_name": "optinist",
    },
    "granger": {
        "function": Granger,
        "conda_name": "optinist",
    },
}
