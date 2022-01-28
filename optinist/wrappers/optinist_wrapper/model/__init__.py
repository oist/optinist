from .glm import GLM
from .cca import CCA
from .svm import SVM
from .lda import LDA


model_wrapper_dict = {
    'glm': {
        'function': GLM
    },
    'cca': {
        'function': CCA
    },
    'svm': {
        'function': SVM
    },
    'lda': {
        'function': LDA
    },
}