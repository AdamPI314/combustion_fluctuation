"""
to deal with numerical interpolation
"""
import time
import numpy as np
import os
import sys


def interp1d(x, y, target):
    """
    return 1d linear interpolation
    """
    return np.interp(target, x, y)


if __name__ == '__main__':
    INIT_TIME = time.time()

    DATA_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print(DATA_DIR)

    x = [1, 2, 3, 4]
    y = [1, 2, 3, 4]
    print(interp1d(x, y, 2.5))

    END_TIME = time.time()

    print("Time Elapsed:\t{:.5} seconds".format(END_TIME - INIT_TIME))
