import numpy as np
from scipy.optimize import bisect

eta = 0.5 / 8.0e3

N = 256

def f( z ):
    return eta * ( z**N - 1.0 ) - ( z - 1.0 )

a = 1.0001
b = 1.8

zoom = bisect( f, a, b, xtol = 1.0e-16 )

print( '{:.15f}'.format( zoom ) )
