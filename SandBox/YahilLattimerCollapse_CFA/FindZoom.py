#!/usr/bin/env python3

import numpy as np
from scipy.optimize import bisect

eta = 0.25 / 1.0e5

N = 512

def f( z ):
    return eta * ( z**N - 1.0 ) - ( z - 1.0 )

a = 1.0001
b = 1.8

zoom = bisect( f, a, b, xtol = 1.0e-16 )

print( '{:.15f}'.format( zoom ) )
