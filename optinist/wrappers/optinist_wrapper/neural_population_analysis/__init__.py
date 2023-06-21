from optinist.wrappers.optinist_wrapper.neural_population_analysis.correlation import (
    correlation,
)
from optinist.wrappers.optinist_wrapper.neural_population_analysis.cross_correlation import (  # noqa: E501
    cross_correlation,
)
from optinist.wrappers.optinist_wrapper.neural_population_analysis.granger import (
    Granger,
)

neural_population_analysis_wrapper_dict = {
    "correlation": {
        "function": correlation,
    },
    "cross_correlation": {
        "function": cross_correlation,
        "conda_yaml": "optinist_env.yaml",
    },
    "granger": {
        "function": Granger,
        "conda_yaml": "optinist_env.yaml",
    },
}
