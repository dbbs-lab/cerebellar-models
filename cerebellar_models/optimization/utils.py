from contextlib import contextmanager
import os
import sys
import matplotlib.pyplot as plt
import numpy as np


@contextmanager
def _suppress_output():
    devnull = open(os.devnull, "w")
    old_out, old_err = sys.stdout, sys.stderr
    sys.stdout, sys.stderr = devnull, devnull
    try:
        yield
    finally:
        sys.stdout, sys.stderr = old_out, old_err
        devnull.close()
