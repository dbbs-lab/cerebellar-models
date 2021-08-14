import os

__version__ = "0.0.1"


class _obj_list(list):
    pass


def templates():
    path = os.path.join(os.path.dirname(__file__), "..")
    return _obj_list([os.path.abspath(path)])
