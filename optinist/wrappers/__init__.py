from .caiman_wrapper import caiman_wrapper_dict
from .suite2p_wrapper import suite2p_wrapper_dict
from .dummy_wrapper import dummy_wrapper_dict

wrapper_dict = {}
wrapper_dict.update(**caiman_wrapper_dict)
wrapper_dict.update(**suite2p_wrapper_dict)
wrapper_dict.update(**dummy_wrapper_dict)