#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def writeData( amrexFile ):

  data = np.loadtxt( amrexFile )

  iiLevel = np.int64( data[:,0] )
  iX1     = np.int64( data[:,1] )
  iX2     = np.int64( data[:,2] )
  iX3     = np.int64( data[:,3] )
  ddX1    = data[:,4]
  dX2     = data[:,5]
  dX3     = data[:,6]
  XX1     = data[:,7]
  X2      = data[:,8]
  X2      = data[:,9]
  ddata   = data[:,10]

  nLevels = np.unique( iiLevel ).shape[0]

  ind    = np.empty( nLevels, object )
  iLevel = np.empty( nLevels, object )
  data   = np.empty( nLevels, object )
  X1     = np.empty( nLevels, object )
  dX1    = np.empty( nLevels, object )

  for i in range( nLevels ):

    ind   [i] = np.where( iiLevel == i )[0]
    iLevel[i] = iiLevel[ind[i]]
    data  [i] = ddata  [ind[i]]
    X1    [i] = XX1    [ind[i]]
    dX1   [i] = ddX1   [ind[i]]

    indd = np.argsort( X1[i] )

    iLevel[i] = np.copy( iLevel[i][indd] )
    data  [i] = np.copy( data  [i][indd] )
    X1    [i] = np.copy( X1    [i][indd] )
    dX1   [i] = np.copy( dX1   [i][indd] )

  return iLevel, X1, dX1, data, ind

amrexFile = 'CF_D_proc000_00000000.dat'
iLevel, X1, dX1, data, ind = writeData( amrexFile )
nLevels = iLevel.shape[0]

for i in range( nLevels-1 ):
  xL = X1[i+1][0:-1]
  xH = X1[i+1][1:  ]
  dX = ( xH - xL )[0]
  for j in range( xL.shape[0] ):
    for k in range( X1[i].shape[0] ):
      if ( ( X1[i][k] > xL[j] ) & ( X1[i][k] < xL[j] + dX ) ) :
        UK0 = data[i][k]
        UK1 = np.sum( data[i+1][j:j+2] ) / 2.0
        AbsValRelDiff = np.abs( UK1 - UK0 ) / ( 0.5 * ( UK1 + UK0 ) )
        print( 'xL, xH, AbsValRelDiff: {:.16e} {:.16e} {:.16e}' \
               .format( xL[j], xH[j], AbsValRelDiff ) )

fig, ax = plt.subplots( 1, 1 )

ms = np.linspace( 1, nLevels+1, nLevels, dtype = int )
for i in range( nLevels ):
  ax.plot( X1[i], data[i], 'ko', ms = ms[i], ls = 'none' )

plt.show()
