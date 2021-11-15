#!/usr/bin/env python3

import yt
import numpy as np
from os import listdir
from os.path import isfile
from sys import exit

yt.funcs.mylog.setLevel(40) # Suppress initial yt output to screen

def MakeDataFile( Field, DataDirectory, DataFileName, TimeFileName, \
                  PlotFileBaseName, UsePhysicalUnits, \
                  WriteExtras = False ):

    if( UsePhysicalUnits ):
        c = 2.99792458e10
    else:
        c = 1.0

    # Get array of all plot-files

    FileArray \
      = np.sort(np.array( [ file for file in listdir( DataDirectory ) ] ) )

    FileList = []

    for iFile in range( FileArray.shape[0] ):

        sFile = FileArray[iFile]

        if( sFile[0:len(PlotFileBaseName)+1] == PlotFileBaseName + '_' \
              and sFile[len(PlotFileBaseName)+1].isdigit() ):
            FileList.append( sFile )

    FileArray = np.array( FileList )

    # Get some general info about the computational domain
    ds       = yt.load( '{:}'.format( DataDirectory + FileArray[0] ) )
    MaxLevel = ds.index.max_level
    nX       = ds.domain_dimensions
    xL       = ds.domain_left_edge
    xH       = ds.domain_right_edge

    if( WriteExtras ):

        with open( 'FileArray.txt', 'w' ) as f:
            for i in FileArray:
                f.write( i )
                f.write( '\n' )
        with open( 'Numbers.txt', 'w' ) as f:
            for i in nX:
                f.write( str(i) )
                f.write( ' ' )
            f.write( '\n' )
            for i in xL.to_ndarray():
                f.write( str(i) )
                f.write( ' ' )
            f.write( '\n' )
            for i in xH.to_ndarray():
                f.write( str(i) )
                f.write( ' ' )
            f.write( '\n' )

        exit()

    Overwrite = True
    if( isfile( DataFileName ) ):

        Overwrite = input( 'File: "{:}" exists. overwrite? (Y/N): '.format \
                      ( DataFileName ) )
        if( not Overwrite == 'Y' ):
            print( 'Not overwriting file' )
            Overwrite = False
        else:
            Overwrite = True

    if( Overwrite ):

        # Put all time-slices into one array to use for movie making
        Data = np.empty( (FileArray.shape[0],nX[0],nX[1]), float )
        Time = np.empty( FileArray.shape[0], float )
        print( 'Generating data file: {:}...'.format( DataFileName ) )
        for i in range( FileArray.shape[0] ):
            print( '{:}/{:}'.format( i+1, FileArray.shape[0] ) )
            ds = yt.load( '{:}'.format( DataDirectory + FileArray[i] ) )

            CoveringGrid \
              = ds.covering_grid \
                  ( level           = MaxLevel, \
                    left_edge       = xL, \
                    dims            = nX * 2**MaxLevel, \
                    num_ghost_zones = nX[0] )

            if( Field == 'PolytropicConstant' ):

                PF_D  \
                  = CoveringGrid['PF_D' ].to_ndarray()[:,:,0]
                AF_P  \
                  = CoveringGrid['AF_P' ].to_ndarray()[:,:,0]
                AF_Gm \
                  = CoveringGrid['AF_Gm'].to_ndarray()[:,:,0]

                Data[i] \
                  = AF_P / PF_D**AF_Gm

            elif( Field == 'LorentzFactor' ):

                PF_V1 \
                  = CoveringGrid['PF_V1'   ].to_ndarray()[:,:,0]
                PF_V2 \
                  = CoveringGrid['PF_V2'   ].to_ndarray()[:,:,0]
                PF_V3 \
                  = CoveringGrid['PF_V3'   ].to_ndarray()[:,:,0]
                GF_g1 \
                  = CoveringGrid['GF_Gm_11'].to_ndarray()[:,:,0]
                GF_g2 \
                  = CoveringGrid['GF_Gm_22'].to_ndarray()[:,:,0]
                GF_g3 \
                  = CoveringGrid['GF_Gm_33'].to_ndarray()[:,:,0]

                Data[i] \
                  = 1 / np.sqrt( 1.0 - ( GF_g1 * PF_V1**2 \
                                         + GF_g2 * PF_V2**2 \
                                         + GF_g3 * PF_V3**2 ) / c**2 )

            elif( Field == 'SpecificEnthalpy' ):

                PF_D \
                  = CoveringGrid['PF_D'].to_ndarray()[:,:,0]
                PF_E \
                  = CoveringGrid['PF_E'].to_ndarray()[:,:,0]
                AF_P \
                  = CoveringGrid['AF_P'].to_ndarray()[:,:,0]

                Data[i] = ( PF_E + AF_P ) / PF_D

            elif( Field == 'Vorticity' ):

                dX1 = ( xH[0].to_ndarray() - xL[0].to_ndarray() ) / nX[0]
                dX2 = ( xH[1].to_ndarray() - xL[1].to_ndarray() ) / nX[1]

                XL = xL.to_ndarray() + 0.5 * np.array( [ dX1, dX2, 0.0 ] )
                XH = xH.to_ndarray() - 0.5 * np.array( [ dX1, dX2, 0.0 ] )

                X1 = np.linspace( XL[0], XH[0], nX[0] )
                X2 = np.linspace( XL[1], XH[1], nX[1] )

                PF_V1 = CoveringGrid['PF_V1'][:,:,0].to_ndarray()
                PF_V2 = CoveringGrid['PF_V2'][:,:,0].to_ndarray()
                indX1 = np.linspace( 1, nX[0]-2, nX[0]-2, dtype = int )
                indX2 = np.linspace( 1, nX[1]-2, nX[1]-2, dtype = int )
                Data[i] = 0.0
                for j in indX1:
                    for k in indX2:
                        Data[i,j,k] \
                          = 1.0 / X1[j] \
                              * ( ( X1[j+1]**2 * PF_V2[j+1,k] \
                                      + X1[j-1]**2 * PF_V2[j-1,k] \
                                      - 2.0 * X1[j]**2 * PF_V2[j,k] ) \
                                  / dX1**2 \
                                    - ( PF_V1[j,k+1] \
                                          + PF_V1[j,k-1] \
                                          - 2.0 * PF_V1[j,k] ) \
                                  / dX2**2 )

            else:

                Data[i] \
                  = CoveringGrid[Field].to_ndarray()[:,:,0]

            Time[i] = ds.current_time

        # Save multi-D array with np.savetxt. Taken from:
        # https://stackoverflow.com/questions/3685265/
        # how-to-write-a-multidimensional-array-to-a-text-file

        with open( DataFileName, 'w' ) as FileOut:
            FileOut.write( '# Array shape: {0}\n'.format( Data.shape ) )

            # Iterating through an n-dimensional array produces slices along
            # the last axis. This is equivalent to Data[i] in this case
            for TimeSlice in Data:
                np.savetxt( FileOut, TimeSlice )
                FileOut.write( '# New slice\n' )

        np.savetxt( TimeFileName, Time )

    return xL, xH, nX, FileArray
