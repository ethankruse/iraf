from .imstat import loadparams
import numpy as np
import matplotlib.pyplot as plt


def implot(*args, **kwargs):
    params = loadparams(*args, **kwargs)

    return params
