from studio.app.optinist.wrappers.optinist.neural_population_analysis.correlation import (  # noqa: E501
    correlation,
)
from studio.app.optinist.wrappers.optinist.neural_population_analysis.cross_correlation import (  # noqa: E501
    cross_correlation,
)
from studio.app.optinist.wrappers.optinist.neural_population_analysis.granger import (
    Granger,
)

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
