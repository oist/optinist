from .cca import CCA
from .pca import PCA
from .tsne import TSNE


dimension_reduction_wrapper_dict = {
    'cca': {
        'function': CCA,
        'conda': 'optinist_env.yaml',
    },
    'pca': {
        'function': PCA,
        'conda': 'optinist_env.yaml',
    },
    'tsne': {
        'function':  TSNE,
        'conda': 'optinist_env.yaml',
    },
}