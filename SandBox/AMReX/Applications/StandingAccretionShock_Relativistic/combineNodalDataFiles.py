#!/usr/bin/env python3

import numpy as np
from datetime import datetime

"""
combineNodalDataFiles.py

Combine files generated by thornado's `ShowVariableFromMultiFab` routine
with multiple processes into single files:
ID_X1n.dat
ID_PF_D.dat

For use in initializing accretion shock problem with multiple levels
"""

nProcs = 4
StepNo = '00000002'
nN     = 3

field = [ 'PF_D', 'PF_V1', 'AF_P' ]

ID = 'GR1D_M2.0_Rpns040_Rsh1.50e2'

### END OF USER INPUT

nLeaves      = np.empty( nProcs, np.int64 )
domainData_r = np.empty( nProcs, object )
domainData_u = np.empty( nProcs, object )

for f in range( len( field ) ):

    nLeafCells = 0

    for i in range( nProcs ):

        data = np.loadtxt( '{:}_{:}_proc{:}_{:}.dat' \
                           .format( ID, field[f], str( i ).zfill( 3 ), StepNo ) )

        ind = np.argsort( data[:,7] )

        rn = np.copy( data[ind,7:7+nN].flatten() )
        un = np.copy( data[ind,-nN:]  .flatten() )

        nLeaves[i] = ind.shape[0]

        nLeafCells += nLeaves[i]

        domainData_r[i] = rn
        domainData_u[i] = un

    combinedData_r = np.empty( nLeafCells*nN, np.float64 )
    combinedData_u = np.empty( nLeafCells*nN, np.float64 )

    iLo = 0
    iHi = 0

    for i in range( nProcs ):

        iHi += nLeaves[i] * nN

        combinedData_r[iLo:iHi] = domainData_r[i]
        combinedData_u[iLo:iHi] = domainData_u[i]

        iLo = iHi

    ind = np.argsort( combinedData_r )

    combinedData_r = np.copy( combinedData_r[ind] )
    combinedData_u = np.copy( combinedData_u[ind] )

    unique = np.unique( combinedData_r, return_index = True )
    ind = unique[1]

    filenameRoot = ID

    if ( field[f] == 'PF_D' ) :
        filename = filenameRoot + '_X1n.dat'
        header \
          = '{:}\nGenerated by {:}\non {:}' \
            .format( filename, __file__, datetime.today() )
        np.savetxt( '{:}'.format( filename ), combinedData_r[ind], \
                    header = header )

    filename = filenameRoot + '_{:}.dat'.format( field[f] )
    header \
      = '{:}\nGenerated by {:}\non {:}' \
        .format( filename, __file__, datetime.today() )

    np.savetxt( '{:}'.format( filename ), combinedData_u[ind], header = header )