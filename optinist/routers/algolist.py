import inspect
from typing import Dict, List, ValuesView

from fastapi import APIRouter

from optinist.routers.const import NOT_DISPLAY_ARGS_LIST
from optinist.routers.model import Algo, AlgoList, Arg, Return
from optinist.wrappers import wrapper_dict

router = APIRouter()


class NestDictGetter:
    @classmethod
    def get_nest_dict(cls, parent_value, parent_key: str) -> Dict[str, Algo]:
        algo_dict = {}
        for key, value in parent_value.items():
            algo_dict[key] = {}
            if isinstance(value, dict) and "function" not in value:
                algo_dict[key]["children"] = cls.get_nest_dict(
                    value, cls._parent_key(parent_key, key)
                )
            else:
                sig = inspect.signature(value["function"])
                returns_list = None
                if sig.return_annotation is not inspect._empty:
                    returns_list = cls._return_list(sig.return_annotation.items())

                algo_dict[key] = Algo(
                    args=cls._args_list(sig.parameters.values()),
                    returns=returns_list,
                    parameter=value["parameter"] if "parameter" in value else None,
                    path=cls._parent_key(parent_key, key),
                )

        return algo_dict

    @classmethod
    def _args_list(cls, arg_params: ValuesView[inspect.Parameter]) -> List[Arg]:
        return [
            Arg(
                name=x.name,
                type=x.annotation.__name__,
                isNone=x.default is None,
            )
            for x in arg_params
            if x.name not in NOT_DISPLAY_ARGS_LIST
        ]

    @classmethod
    def _return_list(cls, return_params: ValuesView[inspect.Parameter]) -> List[Return]:
        return [Return(name=k, type=v.__name__) for k, v in return_params]

    @classmethod
    def _parent_key(cls, parent_key: str, key: str) -> str:
        if parent_key == "":
            return key
        else:
            return f"{parent_key}/{key}"


@router.get("/algolist", response_model=AlgoList, tags=["others"])
async def get_algolist() -> Dict[str, Algo]:
    """_summary_

    Returns:
        {
            'caiman': {
                'children': {
                    'caiman_mc' : {
                        'args': ['images', 'timeseries'],
                        'return': ['images'],
                        'path': 'caiman/caiman_mc'
                    },
                    'caiman_cnmf': {
                        'args': ['images', 'timeseries'],
                        'return': ['images'],
                        'path': 'caiman/caiman_mc'
                    }
                }
            }
        }
    """

    return NestDictGetter.get_nest_dict(wrapper_dict, "")
