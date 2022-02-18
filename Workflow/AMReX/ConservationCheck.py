#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

nN = np.array( [ '02', '03' ], str )
nX = np.array( [ '032' ], str )

Mass = np.ones( (3,nN.shape[0],nX.shape[0]), np.float64 )

Root = '/home/kkadoogan/Work/Codes/thornado/'
Root += 'SandBox/AMReX/Euler_Relativistic_IDEAL/'

def GetData( nN, nX, Grid, Reflux = '' ):

    i = np.loadtxt( Root + 'nN{:}_nX{:}_{:}_I{:}.dat'.format \
          ( nN, nX, Grid, Reflux ) )
    f = np.loadtxt( Root + 'nN{:}_nX{:}_{:}_F{:}.dat'.format \
          ( nN, nX, Grid, Reflux ) )

    N = np.int64( nN )

    iLevel = i[:,0]
    iX1    = i[:,1]
    iX2    = i[:,2]
    iX3    = i[:,3]
    dx1i   = i[:,4]
    dx2    = i[:,5]
    dx3    = i[:,6]
    x1i    = i[:,7:7+N].flatten()
    x2i    = i[:,7+N].flatten()
    x3i    = i[:,7+N+1].flatten()
    ui     = i[:,7+N+2:].flatten()
    di     = np.vstack( (x1i,ui) )
    indi   = di[0].argsort() # sort by x-value

    iLevel = f[:,0]
    iX1    = f[:,1]
    iX2    = f[:,2]
    iX3    = f[:,3]
    dx1f   = f[:,4]
    dx2    = f[:,5]
    dx3    = f[:,6]
    x1f    = f[:,7:7+N].flatten()
    x2f    = f[:,7+N].flatten()
    x3f    = f[:,7+N+1].flatten()
    uf     = f[:,7+N+2:].flatten()
    df     = np.vstack( (x1f,uf) )
    indf   = df[0].argsort() # sort by x-value

    return dx1i, dx1f, di, df, indi, indf

def ComputeMass( N, u, dx ):

    if   N == 1:
        wq = np.array( [ 1.0 ], np.float64 )
    elif N == 2:
        wq = np.array( [ 0.5, 0.5 ], np.float64 )
    elif N == 3:
        wq = np.array( [ 5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0 ], np.float64 )

    Mass = 0.0

    for i in range( dx.shape[0] ):
        jLo = i*N
        jHi = jLo+N
        Mass += np.sum( u[jLo:jHi] * wq ) * dx[i]

    return Mass

for i in range( nN.shape[0] ):

    for j in range( nX.shape[0] ):

#        dxi, dxf, di, df, indi, indf \
#          = GetData( nN[i], nX[j], 'UniGrid' )
#        Mi = ComputeMass( np.int64( nN[i] ), di[1], dxi )
#        Mf = ComputeMass( np.int64( nN[i] ), df[1], dxf )
#        Mass[0,i,j] = abs( ( Mf - Mi ) / Mi )

        dxi, dxf, di, df, indi, indf \
          = GetData( nN[i], nX[j], 'MultiGrid', Reflux = '_RefluxOff' )
        Mi = ComputeMass( np.int64( nN[i] ), di[1], dxi )
        Mf = ComputeMass( np.int64( nN[i] ), df[1], dxf )
        Mass[1,i,j] = abs( ( Mf - Mi ) / Mi )

        dxi, dxf, di, df, indi, indf \
          = GetData( nN[i], nX[j], 'MultiGrid', Reflux = '_RefluxOn' )
        Mi = ComputeMass( np.int64( nN[i] ), di[1], dxi )
        Mf = ComputeMass( np.int64( nN[i] ), df[1], dxf )
        Mass[2,i,j] = abs( ( Mf - Mi ) / Mi )

        for k in range( Mass.shape[0] ):
            Mass[k,i,j] = max( Mass[k,i,j], 1.0e-17 )

nX = np.int64( nX )
nN = np.int64( nN )

fig, ax = plt.subplots( 1, 1 )

for N in range( nN.shape[0] ):
    for K in range( nX.shape[0] ):
        nDOFX = nN[N] * nX[K]
        ax.plot( nDOFX, Mass[0,N,K], 'rs', markersize = 10, \
                 label = 'nN{:}_nX{:}_UniGrid'.format( nN[N], nX[K] ) )
        ax.plot( nDOFX, Mass[1,N,K], 'ms', markersize = 5, \
                 label = 'nN{:}_nX{:}_MultiGrid_RfxOff'.format( nN[N], nX[K] ) )
        ax.plot( nDOFX, Mass[2,N,K], 'bs', markersize = 2, \
                 label = 'nN{:}_nX{:}_MultiGrid_RfxOn'.format( nN[N], nX[K] ) )

ax.legend()
ax.set_xlabel( 'nDOF' )
ax.set_ylabel( r'$\Delta M$' )
ax.set_yscale( 'log' )
plt.show()
#plt.savefig( 'fig.ConservationCheck.png', dpi = 300 )
