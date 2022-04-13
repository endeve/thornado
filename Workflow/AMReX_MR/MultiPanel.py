#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

nN = np.array( [ '01', '02', '03' ], np.str )
nX = np.array( [ '032', '128' ], np.str )

nRows = nN.shape[0]
nCols = nX.shape[0]
fig, axs = plt.subplots( nRows, nCols, figsize = (12,9) )

PlotError = True

if not PlotError:

    for i in range( nRows ):

        for j in range( nCols ):

            X10, dX10, Data0, X11, dX11, Data1 \
              = np.loadtxt( 'nN{:}_nX{:}.dat'.format( nN[i], nX[j] ) )

            if i == 0 and j == 0:
                axs[i,j].plot( X10, Data0, 'r.', label = 't = 0' )
                axs[i,j].plot( X11, Data1, 'k.', label = 't = 10' )
            else:
                axs[i,j].plot( X10, Data0, 'r.' )
                axs[i,j].plot( X11, Data1, 'k.' )
            axs[i,j].set_ylabel( 'nN{:}_nX{:}'.format( nN[i], nX[j] ) )

    axs[0,0].legend()

else:

    fig.suptitle( 'Relative Error' )

    for i in range( nRows ):

        for j in range( nCols ):

            X10U, dX10U, Data0U, X11U, dX11U, Data1U \
              = np.loadtxt( 'nN{:}_nX{:}_UniGrid.dat'.format( nN[i], nX[j] ) )

            X10M, dX10M, Data0M, X11M, dX11M, Data1M \
              = np.loadtxt( 'nN{:}_nX{:}_MultiGrid.dat'.format( nN[i], nX[j] ) )

            RelativeError_Uni   = ( Data1U - Data0U ) / Data0U
            RelativeError_Multi = ( Data1M - Data0M ) / Data0M

            axs[i,j].plot( X10U, RelativeError_Uni, 'r.', label='UniGrid' )
            axs[i,j].plot( X10M, RelativeError_Multi, 'k.', label='MultiGrid' )
            axs[i,j].set_ylabel( 'nN{:}_nX{:}'.format( nN[i], nX[j] ) )

axs[0,0].legend()

plt.subplots_adjust( wspace = 0.5 )
#plt.savefig( 'fig.MultiLevel_RelativeError.png', dpi = 300 )
plt.show()
