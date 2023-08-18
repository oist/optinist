from studio.app.optinist.wrappers.caiman import caiman_wrapper_dict
from studio.app.optinist.wrappers.lccd import lccd_wrapper_dict
from studio.app.optinist.wrappers.optinist import optinist_wrapper_dict
from studio.app.optinist.wrappers.suite2p import suite2p_wrapper_dict

wrapper_dict = {}
wrapper_dict.update(**caiman_wrapper_dict)
wrapper_dict.update(**suite2p_wrapper_dict)
wrapper_dict.update(**lccd_wrapper_dict)
wrapper_dict.update(**optinist_wrapper_dict)
