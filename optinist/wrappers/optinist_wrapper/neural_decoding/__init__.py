from .glm import GLM
from .lda import LDA
from .svm import SVM


neural_decoding_wrapper_dict = {
    'glm': {
        'function': GLM
    },
    'lda': {
        'function': LDA
    },
    'svm': {
        'function': SVM
    },
}