class AlgorithmException(Exception):
    """アルゴリズムの実行時に起きる例外の基底クラス"""

    def __init__(self, message: str):
        super().__init__()
        self._message = message

    def get_message(self) -> str:
        return self._message


class ArgsMissingException(AlgorithmException):
    pass


class ArgsTypeException(AlgorithmException):
    pass
