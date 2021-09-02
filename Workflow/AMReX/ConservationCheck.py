#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

nN = np.array( [ '01', '02', '03' ], np.str )
nX = np.array( [ '032', '128' ], np.str )

MassU = np.empty( (nN.shape[0],nX.shape[0]), np.float64 )
MassM = np.empty( (nN.shape[0],nX.shape[0]), np.float64 )

for i in range( nN.shape[0] ):

    for j in range( nX.shape[0] ):

        X10U, dX10U, Data0U, X11U, dX11U, Data1U \
          = np.loadtxt( 'nN{:}_nX{:}_UniGrid.dat'.format( nN[i], nX[j] ) )

        X10M, dX10M, Data0M, X11M, dX11M, Data1M \
          = np.loadtxt( 'nN{:}_nX{:}_MultiGrid.dat'.format( nN[i], nX[j] ) )

        Mu0 = np.sum( dX10U * Data0U )
        Mu1 = np.sum( dX11U * Data1U )
        MassU[i,j] = abs( Mu1 - Mu0 )

        Mm0 = np.sum( dX10M * Data0M )
        Mm1 = np.sum( dX11M * Data1M )
        MassM[i,j] = abs( Mm1 - Mm0 )

#nX = np.float64( nX )
nN = np.int64( nN )
fig, ax = plt.subplots( 1, 1 )

ax.plot( nN, MassM[:,0], 'k.', label = 'nX = ' + nX[0] + ', MultiGrid' )
ax.plot( nN, MassM[:,1], 'r.', label = 'nX = ' + nX[1] + ', MultiGrid' )
ax.plot( nN, MassU[:,0], 'k^', label = 'nX = ' + nX[0] + ', UniGrid' )
ax.plot( nN, MassU[:,1], 'r^', label = 'nX = ' + nX[1] + ', UniGrid' )

ax.legend()
ax.set_xlabel( 'nN' )
ax.set_ylabel( r'$\Delta M$' )
ax.set_yscale( 'log' )
#plt.show()
plt.savefig( 'fig.ConservationCheck.png', dpi = 300 )
