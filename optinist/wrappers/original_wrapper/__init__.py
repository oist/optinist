from .correlation import correlation
from .pca import PCA
from .granger import Granger
from .glm import GLM

original_wrapper_dict = {
    'original': {
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