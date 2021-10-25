import inspect
import functools


def args_check(func):
    @functools.wraps(func)
    def args_type_check_wrapper(*args, **kwargs):
        """
        引数についてアノテーションで指定した型と一致しているかのチェックを行うデコレータ。
        """
        sig = inspect.signature(func)
        for arg_key, arg_val in sig.bind(*args, **kwargs).arguments.items():
            annot = sig.parameters[arg_key].annotation
            request_type = annot if type(annot) is type else inspect._empty
            if request_type is not inspect._empty and type(arg_val) is not request_type:
                error_msg = f'args"{arg_key}" (type: {type(arg_val)}) is invalid. （ type: {request_type} is require'
                raise TypeError(error_msg)

        return func(*args, **kwargs)

    return args_type_check_wrapper
