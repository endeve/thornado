#!/usr/bin/env python3

import numpy as np

def GetData( dataDirectory, fileNamePrefix, nN ):

    """
    dataDirectory: directory holding the nodal data for each process;
                   e.g., /home/kkadoogan/thornado/SandBox/AMReX/Applications/DynamicTOV_XCFC/DynamicTOV.plt00000000_nodal/
    fileNamePrefix: prefix for nodal data files; e.g., `CF_D` in `CF_D_proc000.dat`
    """

    nProcs = GetNumberOfProcesses( dataDirectory, fileNamePrefix )

    fileNameRoot = dataDirectory + fileNamePrefix

    # Assumes you're working with 1D data
    x1n_G = []
    un_G  = []

    for i in range( nProcs ):

        data \
          = np.loadtxt \
              ( fileNameRoot \
                  + '_proc{:}.dat'.format( str( i ).zfill( 3 ) ) )

        j = 0
        iLevel = np.int64  ( data[:,j]         ); j += 1
        iX1    = np.int64  ( data[:,j]         ); j += 1
        iX2    = np.int64  ( data[:,j]         ); j += 1
        iX3    = np.int64  ( data[:,j]         ); j += 1
        dx1    = np.float64( data[:,j]         ); j += 1
        dx2    = np.float64( data[:,j]         ); j += 1
        dx3    = np.float64( data[:,j]         ); j += 1
        x1n    = np.float64( data[:,j:j+nN[0]] ); j += nN[0]
        x2n    = np.float64( data[:,j:j+nN[1]] ); j += nN[1]
        x3n    = np.float64( data[:,j:j+nN[2]] ); j += nN[2]
        un     = data[:,j:]

        x1n = np.copy( x1n.flatten() )
        x2n = np.copy( x2n.flatten() )
        x3n = np.copy( x3n.flatten() )
        un  = np.copy( un .flatten() )

        x1n_G.append( x1n )
        un_G .append( un )

    flat_x1 = []
    for xProc in x1n_G:
        for x1n in xProc:
            flat_x1.append( x1n )

    flat_u = []
    for uProc in un_G:
        for un in uProc:
            flat_u.append( un )

    # Assumes you're working with 1D data

    x1 = np.array( flat_x1 )
    u  = np.array( flat_u  )

    indX1 = np.argsort( x1 )

    x1n = np.copy( x1[indX1] )
    un  = np.copy( u [indX1] )

    return x1n, un
# END GetData

def GetNumberOfProcesses( dataDirectory, fileNamePrefix ):

    from os import listdir

    allFiles = listdir( dataDirectory )
    nProcs = 0
    for file in allFiles:
        if ( fileNamePrefix in file ) : nProcs += 1

    return nProcs
# END GetNumberOfProcesses


def GetQuadratureWeights1D( nN ):

    if   nN == 1:
        wq = np.array( [ 1.0 ], np.float64 )
    elif nN == 2:
        wq = np.array( [ 0.5, 0.5 ], np.float64 )
    elif nN == 3:
        wq = np.array( [ 5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0 ], np.float64 )

    return wq
# GetQuadratureWeights1D


if ( __name__ == '__main__' ) :

    # User-defined data

    rootDir = '/home/kkadoogan/Work/Codes/thornado/'
    rootDir += 'SandBox/AMReX/Applications/DynamicTOV_XCFC/'

    # Number of DG nodes per element and per dimension
    nN = np.array( [ 2, 1, 1 ], dtype = np.int64 )

    ID     = 'DynamicTOV_nN02_nX032' # PlotFileNameRoot
    StepNo = 101

    fileNamePrefix = 'CF_D'

    # END User-defined data

    dataDirectory \
      = rootDir + '{:}.plt{:}_nodal/'.format( ID, str( StepNo ).zfill( 8 ) )

    x1n, un = GetData( dataDirectory, fileNamePrefix, nN )

    # Unit conversions

    SpeedOfLightMKS          = 2.99792458e8
    GravitationalConstantMKS = 6.673e-11

    Meter    = 1.0
    Second   = SpeedOfLightMKS * Meter
    Kilogram = GravitationalConstantMKS * Meter**3 / Second**2

    Kilometer  = 1.0e3  * Meter
    Centimeter = 1.0e-2 * Meter

    Gram = 1.0e-3 * Kilogram

    MassDensityUnit = Gram / Centimeter**3

    # Plot nodal data

    import matplotlib.pyplot as plt

    plt.plot( x1n / Kilometer, un / MassDensityUnit )
    plt.show()
