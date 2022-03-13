#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from os.path import isfile

Root = '/home/kkadoogan/Work/Codes/thornado/'
Root += 'SandBox/AMReX/Euler_Relativistic_IDEAL/'

Grid = np.array( [ 'MultiGrid' ], str )
nN  = np.array( [ '02' ], str )
nX1 = np.array( [ '016' ], str )
nX2 = np.array( [ '001' ], str )
nX3 = np.array( [ '001' ], str )

shape = (Grid.shape[0],nN.shape[0],nX1.shape[0],nX2.shape[0],nX3.shape[0])

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
                        FileNamesI[g,n,k1,k2,k3] \
                          = Root + 'nN{:}_nX{:}_{:}_I.dat'.format \
                                  ( nN[n], nX1[k1], Grid[g] )
                        FileNamesF[g,n,k1,k2,k3] \
                          = Root + 'nN{:}_nX{:}_{:}_F.dat'.format \
                                  ( nN[n], nX1[k1], Grid[g] )
                    elif( nDims == 2 ):
                        FileNamesI[g,n,k1,k2,k3] \
                          = Root + 'nN{:}_nX{:}x{:}_{:}_I.dat'.format \
                                  ( nN[n], nX1[k1], nX2[k2], Grid[g] )
                        FileNamesF[g,n,k1,k2,k3] \
                          = Root + 'nN{:}_nX{:}x{:}_{:}_F.dat'.format \
                                  ( nN[n], nX1[k1], nX2[k2], Grid[g] )
                    elif( nDims == 3 ):
                        FileNamesI[g,n,k1,k2,k3] \
                          = Root + 'nN{:}_nX{:}x{:}x{:}_{:}_I.dat'.format \
                                  ( nN[n], nX1[k1], nX2[k2], nX3[k3], Grid[g] )
                        FileNamesF[g,n,k1,k2,k3] \
                          = Root + 'nN{:}_nX{:}x{:}x{:}_{:}_F.dat'.format \
                                  ( nN[n], nX1[k1], nX2[k2], nX3[k3], Grid[g] )

def GetData( g, n, k1, k2, k3 ):

    i = np.loadtxt( FileNamesI[g,n,k1,k2,k3] )
    f = np.loadtxt( FileNamesF[g,n,k1,k2,k3] )

    N = np.int64( nN[n] )
    N = np.array( [ 1, 1, 1 ], np.int64 )
    N[0] = nN[n]
    if nDims > 1: N[1] = nN[n]
    if nDims > 2: N[2] = nN[n]

    iLevel = i[:,0]
    iX1    = i[:,1]
    iX2    = i[:,2]
    iX3    = i[:,3]
    dx1i   = i[:,4]
    dx2i   = i[:,5]
    dx3i   = i[:,6]
    x1i    = i[:,7:7+N[0]]
    x2i    = i[:,7+N[0]:7+N[0]+N[1]]
    x3i    = i[:,7+N[0]+N[1]:7+N[0]+N[1]+N[2]]
    ui     = i[:,7+N[0]+N[1]+N[2]:]

    iLevel = f[:,0]
    iX1    = f[:,1]
    iX2    = f[:,2]
    iX3    = f[:,3]
    dx1f   = f[:,4]
    dx2f   = f[:,5]
    dx3f   = f[:,6]
    x1f    = f[:,7:7+N[0]]
    x2f    = f[:,7+N[0]:7+N[0]+N[1]]
    x3f    = f[:,7+N[0]+N[1]:7+N[0]+N[1]+N[2]]
    uf     = f[:,7+N[0]+N[1]+N[2]:]

    return x1i, x2i, x3i, dx1i, dx2i, dx3i, ui, \
           x1f, x2f, x3f, dx1f, dx2f, dx3f, uf

def GetQuadrature1D( nN ):

    if   nN == 1:
        wq = np.array( [ 1.0 ], np.float64 )
    elif nN == 2:
        wq = np.array( [ 0.5, 0.5 ], np.float64 )
    elif N == 3:
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

def ComputeMass( nN, u, dx1, dx2, dx3 ):

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
         Mass += np.sum( u[i] * wq ) * dx1[i] * dx2[i] * dx3[i]

    return Mass

def ComputeL1Error( nN, uf, ui, dx1, dx2, dx3 ):

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

    du = uf - ui

    N = ui.shape[0]
    for i in range( N ):
        L1 += np.sum( np.abs( du[i] ) * wq ) * dx1[i] * dx2[i] * dx3[i]

    return L1

for g in range( Grid.shape[0] ):
    for n in range( nN.shape[0] ):
        for k1 in range( nX1.shape[0] ):
            for k2 in range( nX2.shape[0] ):
                for k3 in range( nX3.shape[0] ):

                    if not ( isfile( FileNamesI[g,n,k1,k2,k3] ) & \
                             isfile( FileNamesF[g,n,k1,k2,k3] ) ): continue

                    x1i, x2i, x3i, dx1i, dx2i, dx3i, ui, \
                    x1f, x2f, x3f, dx1f, dx2f, dx3f, uf \
                      = GetData( g, n, k1, k2, k3 )

                    MassI \
                      = ComputeMass( np.int64( nN[n] ), ui, dx1i, dx2i, dx3i )
                    MassF \
                      = ComputeMass( np.int64( nN[n] ), uf, dx1f, dx2f, dx3f )
                    L1 \
                      = ComputeL1Error \
                          ( np.int64( nN[n] ), uf, ui, dx1i, dx2i, dx3i )
                    print( '\n{:}, N = {:}, nX1 = {:}, nX2 = {:}, nX3 = {:}'.format \
                           ( Grid[g], nN[n], nX1[k1], nX2[k2], nX3[k3] ) )
                    print( '------------------------------------------------' )
                    print( 'dM/M: {:.3e}'.format \
                      ( ( MassF - MassI ) / MassI ) )
                    nDOF = np.int64( np.int64( nN[n] )**nDims )
                    print( 'L1:   {:.3e}'.format \
                      ( L1 / ( nDOF * np.int64( nX1[k1] ) \
                                 * np.int64( nX2[k2] ) * np.int64( nX3[k3] ) ) ) )
