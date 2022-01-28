from .cca import CCA
from .pca import PCA
from .tsne import TSNE


dimension_reduction_wrapper_dict = {
    'cca': {
        'function': CCA
    },
    'pca': {
        'function': PCA
    },
    'tsne': {
        'function':  TSNE
    },
}