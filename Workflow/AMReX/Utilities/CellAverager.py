#!/usr/bin/env python3

import numpy as np

import GlobalVariables.Settings as gvS
import GlobalVariables.Units    as gvU



# Pack all field data into NodalData.
# Use length of NodalData to establish loop over fields.
# Calc cell averages of data.
# Update X1

def computeCellAverage_Old( Names, nNodes, Field, nSS ):

  x1_C = np.empty( nX1, np.float64 )

  uh = Names[Field][1]

  if nNodes == 2:
    wq = np.array( [ 0.5, 0.5 ], np.float64 )
  else:
    exit( 'Not available for nNodes = {:}'.format( nNodes ) )

  SqrtGm = Names['GF_Sg'][1]

  uK = np.empty( (nSS,nX1), np.float64 )

  for iSS in range( nSS ):
    for iX1 in range( nX1 ):

      iLo = nNodes * iX1
      iHi = iLo + nNodes

      vK = np.sum( wq * SqrtGm[iSS,0,0,iLo:iHi] )

      uK[iSS,iX1] \
        = np.sum( wq * uh[iSS,0,0,iLo:iHi] * SqrtGm[iSS,0,0,iLo:iHi] ) / vK

      if iSS == 0:
        x1_C[iX1] = np.sum( wq * x1q[iLo:iHi] )

  return uK, x1_C



def computeCellAverage( NodalData, SqrtGm, NodalX1, nX  ):
    
    nNodes = int(len(NodalData)/nX)
    
    if nNodes == 2:
        wq = np.array( [ 0.5, 0.5 ], np.float64 )
    else:
        exit( 'Not available for nNodes = {:}'.format( nNodes ) )
    
    AverageData = np.empty( (nX), np.float64 )
    AverageX1 = np.empty( (nX), np.float64 )
    for iX in range( nX ):
        
        iLo = nNodes*iX
        iHi = iLo + nNodes
        
        vK = np.sum( wq * SqrtGm[iLo:iHi] )
        
        AverageData[iX] \
        = np.sum( wq * NodalData[iLo:iHi] * SqrtGm[iLo:iHi] ) / vK
        
        AverageX1[iX] \
        = np.sum( wq * NodalX1[iLo:iHi])
    
    
#    exit('Here')
    return AverageData, AverageX1

