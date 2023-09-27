from studio.app.optinist.wrappers.optinist.neural_decoding.glm import GLM
from studio.app.optinist.wrappers.optinist.neural_decoding.lda import LDA
from studio.app.optinist.wrappers.optinist.neural_decoding.svm import SVM

neural_decoding_wrapper_dict = {
    "glm": {
        "function": GLM,
        "conda_name": "optinist",
    },
    "lda": {
        "function": LDA,
        "conda_name": "optinist",
    },
    "svm": {
        "function": SVM,
        "conda_name": "optinist",
    },
}
