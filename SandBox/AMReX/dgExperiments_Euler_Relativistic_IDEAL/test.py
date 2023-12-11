#!/usr/bin/env python3

import numpy as np

from gaussxw import gaussxw

N = 2

xq, wq = gaussxw( N )

def Lagrange( x, i, xx ):
  L = 1.0
  for j in range( xx.shape[0] ):
    if j != i:
      L *= ( x - xx[j] ) / ( xx[i] - xx[j] )
  return L

def r( eta, rC, dr ):
    return rC + eta * dr

rC = np.array( [ 0.5, 1.5 ] )
dr = 1.0
RC = 1.0
dR = 2.0

SqrtGmj = np.array( [ r( xq, rC[0], dr )**2, r( xq, rC[1], dr )**2 ] )
SqrtGm  = r( xq, RC, dR )
print(SqrtGmj)
print(SqrtGm)

