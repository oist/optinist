from .correlation import correlation
from .granger import Granger
from .pca import PCA
from .tsne import TSNE


stats_wrapper_dict = {
    'correlation': {
        'function': correlation
    },
    'granger': {
        'function': Granger
    },
    'pca': {
        'function': PCA
    },
    'tsne': {
        'function':  TSNE
    },
}