from optinist.wrappers.optinist_wrapper.neural_decoding.glm import GLM
from optinist.wrappers.optinist_wrapper.neural_decoding.lda import LDA
from optinist.wrappers.optinist_wrapper.neural_decoding.svm import SVM

neural_decoding_wrapper_dict = {
    "glm": {
        "function": GLM,
        "conda_yaml": "optinist_env.yaml",
    },
    "lda": {
        "function": LDA,
        "conda_yaml": "optinist_env.yaml",
    },
    "svm": {
        "function": SVM,
        "conda_yaml": "optinist_env.yaml",
    },
}
