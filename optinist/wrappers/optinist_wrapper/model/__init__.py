from .glm import GLM
from .cca import CCA

model_wrapper_dict = {
    'glm': {
        'function': GLM
    },
    'cca': {
        'function':  CCA
    },
}