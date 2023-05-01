from optinist.wrappers.optinist_wrapper.dimension_reduction.cca import CCA
from optinist.wrappers.optinist_wrapper.dimension_reduction.pca import PCA
from optinist.wrappers.optinist_wrapper.dimension_reduction.tsne import TSNE

dimension_reduction_wrapper_dict = {
    "cca": {
        "function": CCA,
        "conda_yaml": "optinist_env.yaml",
    },
    "pca": {
        "function": PCA,
        "conda_yaml": "optinist_env.yaml",
    },
    "tsne": {
        "function": TSNE,
        "conda_yaml": "optinist_env.yaml",
    },
}
