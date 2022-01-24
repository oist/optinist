import inspect
import functools

from pynwb import NWBFile
from wrappers.optinist_exception import ArgsMissingException, ArgsTypeException

def args_check(func):
    @functools.wraps(func)
    def args_type_check_wrapper(*args, **kwargs):
        sig = inspect.signature(func)
        args = list(args)

        arg_type_list = [type(x) for x in args]
        request_type_list = [x.annotation for x in sig.parameters.values()]
        request_default_list = [x.default for x in sig.parameters.values()]

        # 引数が足らない場合はエラーさせる
        for x, de in zip(request_type_list, request_default_list):
            if de is not None and arg_type_list.count(x) < request_type_list.count(x):
                error_msg = f'args type error: need " {x.__name__} " type, \
                            more than " {request_type_list.count(x) - arg_type_list.count(x)} " '
                raise ArgsMissingException(error_msg)

        # いらない型の引数を削除する
        for i, x in reversed(list(enumerate(arg_type_list))):
            if x not in request_type_list:
                args.pop(i)

        arg_type_list = [type(x) for x in args]
        request_type_list = [x.annotation for x in sig.parameters.values()]
        request_default_list = [x.default for x in sig.parameters.values()]

        # 型の順番を揃える
        def check_order(arg_type_list, request_type_list, request_default_list):
            not_good = False
            for i in range(len(request_type_list)):
                if len(arg_type_list) > i and arg_type_list[i] != request_type_list[i] and request_type_list[i] != dict:
                    not_good = True
                    break

            return not_good

        count_loop = 0
        while check_order(arg_type_list, request_type_list, request_default_list):
            count_loop += 1
            if count_loop > 100:
                raise "loop error"
            for i in range(len(args)-1):
                if arg_type_list[i] == request_type_list[i]:
                    # import pdb; pdb.set_trace()
                    continue

                if check_order(arg_type_list, request_type_list, request_default_list):
                    args[i], args[i+1] = args[i+1], args[i]
                    arg_type_list[i], arg_type_list[i+1] = arg_type_list[i+1], arg_type_list[i]

            arg_type_list = [type(x) for x in args]
            request_type_list = [x.annotation for x in sig.parameters.values()]

        # 引数が多い場合は余分な引数を捨てる
        args = args[:len(sig.parameters)]
        arg_type_list = [type(x) for x in args]

        # paramsを抜く
        if len(arg_type_list) == len(request_type_list) and \
                request_type_list[-1] == dict and \
                arg_type_list[-1] != request_type_list[-1]:
            args = args[:-1]

        arg_type_list = [type(x) for x in args]
        request_type_list = [x.annotation for x in sig.parameters.values()]

        for arg_key, arg_val in sig.bind(*args, **kwargs).arguments.items():
            # 渡した引数
            arg_type = type(arg_val)
            
            # 関数で定義された引数
            annot = sig.parameters[arg_key].annotation
            request_type = annot if type(annot) is type else inspect._empty
            if request_type is not inspect._empty and arg_type is not request_type:
                error_msg = f'args"{arg_key}" (type: {arg_type}) is invalid. （ type: {request_type} is require'
                raise ArgsTypeException(error_msg)

        return func(*args, **kwargs)

    return args_type_check_wrapper

