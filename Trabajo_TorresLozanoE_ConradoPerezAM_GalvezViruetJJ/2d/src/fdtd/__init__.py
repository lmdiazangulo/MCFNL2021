from __future__ import division

import math
import numpy as np

def isClose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def isLower(x, y):
    if isClose(x,y):
        return False
    else:
        return x < y

def isGreater(x, y):
    if isClose(x,y):
        return False
    else:
        return x > y
