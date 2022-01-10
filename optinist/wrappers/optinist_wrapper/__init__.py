from .correlation import correlation
from .pca import PCA
from .granger import Granger
from .glm import GLM

optinist_wrapper_dict = {
    'optinist': {
        'correlation': {
			'function': correlation
        },
        'pca': {
			'function': PCA
        },
        'granger': {
            'function': Granger
        },
        'glm': {
            'function': GLM
        }
    }
}