from optinist.wrappers.optinist.dimension_reduction.cca import CCA
from optinist.wrappers.optinist.dimension_reduction.pca import PCA
from optinist.wrappers.optinist.dimension_reduction.tsne import TSNE

dimension_reduction_wrapper_dict = {
    "cca": {
        "function": CCA,
        "conda_name": "optinist",
    },
    "pca": {
        "function": PCA,
        "conda_name": "optinist",
    },
    "tsne": {
        "function": TSNE,
        "conda_name": "optinist",
    },
}
