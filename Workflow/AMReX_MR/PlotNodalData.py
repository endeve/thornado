#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from os.path import isfile

Root = '/home/kkadoogan/Work/Codes/thornado/'
Root += 'SandBox/AMReX/Euler_Relativistic_IDEAL_MR/'

nF = 2 # 0: fluid; 1: SqrtGm
Grid = np.array( [ 'UniGrid', 'MultiGrid' ], str)
nN  = np.array( [ '02', '03' ], str )
nX1 = np.array( [ '032', '128' ], str )
nX2 = np.array( [ '001' ], str )
nX3 = np.array( [ '001' ], str )

shape = (nF,Grid.shape[0],nN.shape[0],nX1.shape[0],nX2.shape[0],nX3.shape[0])

FileNamesI = np.empty( shape, object )
FileNamesF = np.empty( shape, object )

# --- FileName depends on dimensionality ---

nDims = 1

for g in range( Grid.shape[0] ):
    for n in range( nN.shape[0] ):
        for k1 in range( nX1.shape[0] ):
            for k2 in range( nX2.shape[0] ):
                for k3 in range( nX3.shape[0] ):
                    if( nDims == 1 ):
                        FileNamesI[0,g,n,k1,k2,k3] \
                          = Root + 'nN{:}_nX{:}_{:}_U_I.dat'.format \
                                  ( nN[n], nX1[k1], Grid[g] )
                        FileNamesF[0,g,n,k1,k2,k3] \
                          = Root + 'nN{:}_nX{:}_{:}_U_F.dat'.format \
                                  ( nN[n], nX1[k1], Grid[g] )
                        FileNamesI[1,g,n,k1,k2,k3] \
                          = Root + 'nN{:}_nX{:}_{:}_G_I.dat'.format \
                                  ( nN[n], nX1[k1], Grid[g] )
                        FileNamesF[1,g,n,k1,k2,k3] \
                          = Root + 'nN{:}_nX{:}_{:}_G_F.dat'.format \
                                  ( nN[n], nX1[k1], Grid[g] )
                    elif( nDims == 2 ):
                        FileNamesI[0,g,n,k1,k2,k3] \
                          = Root + 'nN{:}_nX{:}x{:}_{:}_U_I.dat'.format \
                                  ( nN[n], nX1[k1], nX2[k2], Grid[g] )
                        FileNamesF[0,g,n,k1,k2,k3] \
                          = Root + 'nN{:}_nX{:}x{:}_{:}_U_F.dat'.format \
                                  ( nN[n], nX1[k1], nX2[k2], Grid[g] )
                        FileNamesI[1,g,n,k1,k2,k3] \
                          = Root + 'nN{:}_nX{:}x{:}_{:}_G_I.dat'.format \
                                  ( nN[n], nX1[k1], nX2[k2], Grid[g] )
                        FileNamesF[1,g,n,k1,k2,k3] \
                          = Root + 'nN{:}_nX{:}x{:}_{:}_G_F.dat'.format \
                                  ( nN[n], nX1[k1], nX2[k2], Grid[g] )
                    elif( nDims == 3 ):
                        FileNamesI[0,g,n,k1,k2,k3] \
                          = Root + 'nN{:}_nX{:}x{:}x{:}_{:}_U_I.dat'.format \
                                  ( nN[n], nX1[k1], nX2[k2], nX3[k3], Grid[g] )
                        FileNamesF[0,g,n,k1,k2,k3] \
                          = Root + 'nN{:}_nX{:}x{:}x{:}_{:}_U_F.dat'.format \
                                  ( nN[n], nX1[k1], nX2[k2], nX3[k3], Grid[g] )
                        FileNamesI[1,g,n,k1,k2,k3] \
                          = Root + 'nN{:}_nX{:}x{:}x{:}_{:}_G_I.dat'.format \
                                  ( nN[n], nX1[k1], nX2[k2], nX3[k3], Grid[g] )
                        FileNamesF[1,g,n,k1,k2,k3] \
                          = Root + 'nN{:}_nX{:}x{:}x{:}_{:}_G_F.dat'.format \
                                  ( nN[n], nX1[k1], nX2[k2], nX3[k3], Grid[g] )

def GetData( g, n, k1, k2, k3 ):

    iU = np.loadtxt( FileNamesI[0,g,n,k1,k2,k3] )
    fU = np.loadtxt( FileNamesF[0,g,n,k1,k2,k3] )
    iG = np.loadtxt( FileNamesI[1,g,n,k1,k2,k3] )
    fG = np.loadtxt( FileNamesF[1,g,n,k1,k2,k3] )

    N = np.int64( nN[n] )
    N = np.array( [ 1, 1, 1 ], np.int64 )
    N[0] = nN[n]
    if nDims > 1: N[1] = nN[n]
    if nDims > 2: N[2] = nN[n]

    iLevel = iU[:,0]
    iX1    = iU[:,1]
    iX2    = iU[:,2]
    iX3    = iU[:,3]
    dx1i   = iU[:,4]
    dx2i   = iU[:,5]
    dx3i   = iU[:,6]
    x1i    = iU[:,7:7+N[0]]
    x2i    = iU[:,7+N[0]:7+N[0]+N[1]]
    x3i    = iU[:,7+N[0]+N[1]:7+N[0]+N[1]+N[2]]
    ui     = iU[:,7+N[0]+N[1]+N[2]:]
    gi     = iG[:,7+N[0]+N[1]+N[2]:]

    iLevel = fU[:,0]
    iX1    = fU[:,1]
    iX2    = fU[:,2]
    iX3    = fU[:,3]
    dx1f   = fU[:,4]
    dx2f   = fU[:,5]
    dx3f   = fU[:,6]
    x1f    = fU[:,7:7+N[0]]
    x2f    = fU[:,7+N[0]:7+N[0]+N[1]]
    x3f    = fU[:,7+N[0]+N[1]:7+N[0]+N[1]+N[2]]
    uf     = fU[:,7+N[0]+N[1]+N[2]:]
    gf     = fG[:,7+N[0]+N[1]+N[2]:]

    return x1i, x2i, x3i, dx1i, dx2i, dx3i, ui, gi, \
           x1f, x2f, x3f, dx1f, dx2f, dx3f, uf, gf

def GetQuadrature1D( nN ):

    if   nN == 1:
        wq = np.array( [ 1.0 ], np.float64 )
    elif nN == 2:
        wq = np.array( [ 0.5, 0.5 ], np.float64 )
    elif nN == 3:
        wq = np.array( [ 5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0 ], np.float64 )

    return wq

def GetQuadrature( nN ):

    if nDims == 1:
        wq1 = GetQuadrature1D( nN )
        wq2 = GetQuadrature1D( 1 )
        wq3 = GetQuadrature1D( 1 )
    elif nDims == 2:
        wq1 = GetQuadrature1D( nN )
        wq2 = GetQuadrature1D( nN )
        wq3 = GetQuadrature1D( 1 )
    elif nDims == 3:
        wq1 = GetQuadrature1D( nN )
        wq2 = GetQuadrature1D( nN )
        wq3 = GetQuadrature1D( nN )

    nDOFX = nN**nDims

    wq = np.empty( nDOFX, np.float64 )

    ii = -1
    for i3 in range( wq3.shape[0] ):
        for i2 in range( wq2.shape[0] ):
            for i1 in range( wq1.shape[0] ):
                ii += 1
                wq[ii] = wq1[i1] * wq2[i2] * wq3[i3]

    return wq

def ComputeMass( nN, u, g, dx1, dx2, dx3 ):

    if nDims == 1:
        wq1 = GetQuadrature1D( nN )
        wq2 = GetQuadrature1D( 1 )
        wq3 = GetQuadrature1D( 1 )
    elif nDims == 2:
        wq1 = GetQuadrature1D( nN )
        wq2 = GetQuadrature1D( nN )
        wq3 = GetQuadrature1D( 1 )
    elif nDims == 3:
        wq1 = GetQuadrature1D( nN )
        wq2 = GetQuadrature1D( nN )
        wq3 = GetQuadrature1D( nN )

    nDOFX = nN**nDims

    wq = GetQuadrature( nN )

    Mass = 0.0

    N = u.shape[0]
    for i in range( N ):
         Mass += np.sum( wq * u[i] * g[i] ) * dx1[i] * dx2[i] * dx3[i]

    return Mass

def ComputeL1Error( nN, uf, ui, gf, gi, dx1, dx2, dx3 ):

    if nDims == 1:
        wq1 = GetQuadrature1D( nN )
        wq2 = GetQuadrature1D( 1 )
        wq3 = GetQuadrature1D( 1 )
    elif nDims == 2:
        wq1 = GetQuadrature1D( nN )
        wq2 = GetQuadrature1D( nN )
        wq3 = GetQuadrature1D( 1 )
    elif nDims == 3:
        wq1 = GetQuadrature1D( nN )
        wq2 = GetQuadrature1D( nN )
        wq3 = GetQuadrature1D( nN )

    nDOFX = nN**nDims

    wq = GetQuadrature( nN )

    L1 = 0.0

    du = uf * gf - ui * gi

    N = ui.shape[0]
    for i in range( N ):
        L1 += np.sum( wq * np.abs( du[i] ) ) * dx1[i] * dx2[i] * dx3[i]

    return L1

for g in range( Grid.shape[0] ):
    for n in range( nN.shape[0] ):
        for k1 in range( nX1.shape[0] ):
            for k2 in range( nX2.shape[0] ):
                for k3 in range( nX3.shape[0] ):

                    if not ( isfile( FileNamesI[0,g,n,k1,k2,k3] ) & \
                             isfile( FileNamesF[0,g,n,k1,k2,k3] ) ): continue

                    x1i, x2i, x3i, dx1i, dx2i, dx3i, ui, gi, \
                    x1f, x2f, x3f, dx1f, dx2f, dx3f, uf, gf \
                      = GetData( g, n, k1, k2, k3 )
                    plt.plot( x1i.flatten(), ui.flatten(), 'r.' )
                    plt.plot( x1f.flatten(), uf.flatten(), 'b.' )
                    plt.show()

                    N = np.int64( nN[n] )

                    MassI \
                      = ComputeMass( N, ui, gi, dx1i, dx2i, dx3i )
                    MassF \
                      = ComputeMass( N, uf, gf, dx1f, dx2f, dx3f )
                    L1 \
                      = ComputeL1Error( N, uf, ui, gf, gi, dx1i, dx2i, dx3i )
                    print( '\n{:}, N = {:}, nX1 = {:}, nX2 = {:}, nX3 = {:}'.format \
                           ( Grid[g], nN[n], nX1[k1], nX2[k2], nX3[k3] ) )
                    print( '------------------------------------------------' )
                    print( 'dM/M: {:.3e}'.format \
                      ( ( MassF - MassI ) / MassI ) )
                    nDOF = np.int64( np.int64( nN[n] )**nDims )
                    print( 'L1:   {:.3e}'.format \
                      ( L1 / ( nDOF * np.int64( nX1[k1] ) \
                                 * np.int64( nX2[k2] ) * np.int64( nX3[k3] ) ) ) )
