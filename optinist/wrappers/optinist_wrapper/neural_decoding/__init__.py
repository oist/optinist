from .glm import GLM
from .lda import LDA
from .svm import SVM


neural_decoding_wrapper_dict = {
    'glm': {
        'function': GLM,
        'conda_yaml': 'optinist_env.yaml',
    },
    'lda': {
        'function': LDA,
        'conda_yaml': 'optinist_env.yaml',
    },
    'svm': {
        'function': SVM,
        'conda_yaml': 'optinist_env.yaml',
    },
}