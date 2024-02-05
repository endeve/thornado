#!/usr/bin/env python3

import numpy as np
import h5py as h5
import yt
import math

import GlobalVariables.Settings as gvS

#=============================================#
#   Included Routines
#
#   CreateLocations
#   CreateLeafLocations
#   CreateUniformLocations
#
#=============================================#



 #=============================================#
#                                               #
#   CreateLocations                             #
#                                               #
 #=============================================#
def CreateLocations(File,               \
                    MaxLevel = -1,      \
                    TypeIn = "AMReX"     ):

    CheckMeshType(TypeIn)

    if TypeIn.lower() == "amrex":
        X1, X2, X3, dX1, dX2, dX3, xL, xH           \
            = CreateLeafLocations(  File            )

    elif TypeIn.lower() == "native":
        X1, X2, X3, dX1, dX2, dX3, xL, xH           \
            = CreateNativeLocations(  File          )

    elif TypeIn.lower() == "uniform":
        X1, X2, X3, dX1, dX2, dX3, xL, xH           \
            = CreateUniformLocations(   File,       \
                                        MaxLevel    )


    return X1, X2, X3, dX1, dX2, dX3, xL, xH






 #=============================================#
#                                               #
#   CreateLeafLocations                         #
#                                               #
 #=============================================#
def CreateLeafLocations( File ):

    yt.funcs.mylog.setLevel(40) # Suppress yt warnings
    ds = yt.load( '{:}'.format( File ) )

    Time = ds.current_time.to_ndarray()
    nX   = ds.domain_dimensions
    xL   = ds.domain_left_edge
    xH   = ds.domain_right_edge

    nDimsX = 1
    if nX[1] > 1: nDimsX += 1
    if nX[2] > 1: nDimsX += 1

    X1  = []
    dX1 = []
    for lvl in range(ds.max_level+1):

        grids = ds.index.select_grids(lvl)

        for grid in grids:
            for cell in range(grid.ActiveDimensions[0]):
                X1_C = grid["boxlib","X1_C"]
                dX1g = grid["boxlib","dX1"]

                if grid.child_mask[cell]:
                    X1 .append( X1_C[cell].value[0][0] )
                    dX1.append( dX1g[cell].value[0][0] )

    if   ( gvS.CoordinateSystem.lower() == 'cartesian' ) :

      X2  = [ 0.5 ]
      X3  = [ 0.5 ]
      dX2 = [ 1.0 ]
      dX3 = [ 1.0 ]

    elif ( gvS.CoordinateSystem.lower() == 'spherical' ) :

      X2  = [ math.pi/2.0 ]
      X3  = [ math.pi ]
      dX2 = [ math.pi ]
      dX3 = [ 2.0*math.pi ]

    X1  = np.array( X1  )
    X2  = np.array( X2  )
    X3  = np.array( X3  )
    dX1 = np.array( dX1 )
    dX2 = np.array( dX2 )
    dX3 = np.array( dX3 )

    indX1 = np.argsort( X1 )
    indX2 = np.argsort( X2 )
    indX3 = np.argsort( X3 )

    X1  = np.copy( X1 [indX1] )
    X2  = np.copy( X2 [indX2] )
    X3  = np.copy( X3 [indX3] )

    dX1 = np.copy( dX1[indX1] )
    dX2 = np.copy( dX2[indX2] )
    dX3 = np.copy( dX3[indX3] )

    return X1, X2, X3, dX1, dX2, dX3, xL, xH





 #=============================================#
#                                               #
#   CreateNativeLocations                       #
#                                               #
 #=============================================#
def CreateNativeLocations( File ):

    PathRoot = File[:-6]
    FrameNumber = File[-6:]
    FF_root = PathRoot + 'FluidFields_'

    DataFileName = FF_root + str( FrameNumber ) + '.h5'
    Data_FF      = h5.File( DataFileName, 'r' )

    # Get the spatial grid
    X    = Data_FF[ 'Spatial Grid' ]
    X1   = np.array( X[ 'X1' ] )
    X2   = np.array( X[ 'X2' ] )
    X3   = np.array( X[ 'X3' ] )
    X1_C = np.array( X[ 'X1_C' ] )
    X2_C = np.array( X[ 'X2_C' ] )
    X3_C = np.array( X[ 'X3_C' ] )
    dX1  = np.array( X[ 'dX1' ] )
    dX2  = np.array( X[ 'dX2' ] )
    dX3  = np.array( X[ 'dX3' ] )

    xL = [X1_C[0]-dX1[0]/2.0,X2_C[0]-dX2[0]/2.0,X3_C[0]-dX3[0]/2.0]
    xH = [X1_C[-1]+dX1[-1]/2.0,X2_C[-1]+dX2[-1]/2.0,X3_C[-1]+dX3[-1]/2.0]

    return X1, X2, X3, dX1, dX2, dX3, xL, xH





 #=============================================#
#                                               #
#   CreateUniformLocations                      #
#                                               #
 #=============================================#
def CreateUniformLocations( File,               \
                            MaxLevel            ):


    yt.funcs.mylog.setLevel(40) # Suppress yt warnings
    ds = yt.load( '{:}'.format( File ) )

    if MaxLevel == -1: MaxLevel = ds.index.max_level
    Time = ds.current_time.to_ndarray()
    nX   = ds.domain_dimensions
    xL   = ds.domain_left_edge
    xH   = ds.domain_right_edge

    nDimsX = 1
    if nX[1] > 1: nDimsX += 1
    if nX[2] > 1: nDimsX += 1

    if nDimsX == 1:
        dim_array = [nX[0] * 2**MaxLevel, nX[1], nX[2]]
    elif nDimsX == 2:
        dim_array = [nX[0] * 2**MaxLevel, nX[1] * 2**MaxLevel, nX[2]]
    else:
        dim_array = nX * 2**MaxLevel

    CoveringGrid \
          = ds.covering_grid \
              ( level           = MaxLevel, \
                left_edge       = xL, \
                dims            = dim_array, \
                num_ghost_zones = 0 )




    X1 = np.copy( CoveringGrid['X1_C'].to_ndarray() )
    X1 = np.copy( X1[:,0,0] )
    dX1 = np.copy( CoveringGrid['dX1'].to_ndarray() )
    dX1 = np.copy( dX1[:,0,0])

    if nDimsX > 1:
        X2 = np.copy( CoveringGrid['X2_C'].to_ndarray() )
        dX2 = np.copy( CoveringGrid['dX2'].to_ndarray() )
        X2 = np.copy( X2[0,:,0] )
        dX2 = np.copy( dX2[0,:,0] )
    else:
        X2  = np.array( [math.pi/2.0] )
        dX2 = np.array( [math.pi] )

    if nDimsX > 2:
        X3 = np.copy( CoveringGrid['X3_C'].to_ndarray() )
        dX3 = np.copy( CoveringGrid['dX3'].to_ndarray() )
        X3 = np.copy( X3[0,0,:] )
        dX3 = np.copy( dX3[0,0,:] )
    else:
        X3  = np.array( [math.pi] )
        dX3 = np.array( [2.0*math.pi] )



    del ds, CoveringGrid

    return X1, X2, X3, dX1, dX2, dX3, xL, xH





 #=============================================#
#                                               #
#   CheckMeshType                               #
#                                               #
 #=============================================#
def CheckMeshType( TypeIn ):

    Types = set(Type.lower() for Type in ("amrex","native","uniform"))

    msg =  "\n\n Invalid mesh type: {:}. \n".format( TypeIn )
    msg += "\n Acceptable types: \n"
    msg += "-------------------\n"
    msg += "AMReX \n"
    msg += "Native \n"
    msg += "Uniform \n"
    assert( TypeIn.lower() in Types ), msg

    return
