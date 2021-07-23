#!/usr/bin/env python3

import yt
import numpy as np

yt.funcs.mylog.setLevel(40) # Suppress initial yt output to screen

from UtilitiesModule import OverwriteFile, GetData, ChoosePlotFile

def MakeDataFile( Field, DataDirectory, DataFileName, \
                  PlotFileBaseName, CoordinateSystem, \
                  UsePhysicalUnits = True, ReturnMesh = False, \
                  WriteExtras = False, Verbose = False ):

    print( '\nRunning MakeDataFile...\n' )

    print( 'DataDirectory: {:}\n'.format( DataDirectory ) )

    if( not DataFileName[0] == '.' ): DataFileName = '.' + DataFileName

    c = 2.99792458e10
    if( not UsePhysicalUnits ): c = 1.0

    # Append "/" to DataDirectory, if not present
    if( not DataDirectory[-1] == '/' ): DataDirectory += '/'

    File, FileArray \
      = ChoosePlotFile( DataDirectory, PlotFileBaseName, \
                        ReturnFileArray = True, Verbose = Verbose )

    # Get some general info about the computational domain

    ds       = yt.load( DataDirectory + FileArray[0] )
    MaxLevel = ds.index.max_level
    nX       = ds.domain_dimensions
    xL       = ds.domain_left_edge
    xU       = ds.domain_right_edge

    nDimsX = 1
    if nX[1] > 1: nDimsX += 1
    if nX[2] > 1: nDimsX += 1

    if nDimsX == 3:

      msg = 'MakeDataFile not implemented for nDimsX = 3'

      assert 0, msg

    if WriteExtras:

        # This is if you're running on a computing cluster and don't want
        # to copy the files to your local machine

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
            for i in xU.to_ndarray():
                f.write( str(i) )
                f.write( ' ' )
            f.write( '\n' )

        exit()

    OW = OverwriteFile( DataFileName )

    if OW:

        # Put all time-slices into one array to use for movie making

        if   nDimsX == 1: DataShape = (FileArray.shape[0],nX[0])
        elif nDimsX == 2: DataShape = (FileArray.shape[0],nX[0],nX[1])

        Data = np.empty( DataShape, np.float64 )
        Time = np.empty( FileArray.shape[0], np.float64 )

        print( 'Generating data file: {:}...'.format( DataFileName ) )

        for i in range( FileArray.shape[0] ):

            print( '{:}/{:}'.format( i+1, FileArray.shape[0] ) )

            ds = yt.load( '{:}'.format( DataDirectory + FileArray[i] ) )

            Data[i], DataUnit, Time[i] \
              = GetData( DataDirectory, PlotFileBaseName, \
                         Field, CoordinateSystem, UsePhysicalUnits, \
                         argv = [ 'a', FileArray[i] ], \
                         ReturnTime = True, ReturnMesh = False, \
                         Verbose = Verbose )

        # Save multi-D array with np.savetxt. Taken from:
        # https://stackoverflow.com/questions/3685265/
        # how-to-write-a-multidimensional-array-to-a-text-file

        TimeHeader = '# Time [ms] '
        if not UsePhysicalUnits: TimeHeader = '# Time '

        with open( DataFileName, 'w' ) as FileOut:

            FileOut.write( '# Array shape: {:}\n'.format( DataShape ) )
            FileOut.write( '# Units: {:}\n'.format( DataUnit ) )

        with open( DataFileName, 'a' ) as FileOut:

            FileOut.write( TimeHeader )
            for t in Time:
                FileOut.write( str( t ) + ' ' )
            FileOut.write( '\n' )

        with open( DataFileName, 'a' ) as FileOut:

            # Iterating through an n-dimensional array produces slices along
            # the last axis. This is equivalent to Data[i] in this case

            for TimeSlice in Data:
                FileOut.write( '# New slice\n' )
                np.savetxt( FileOut, TimeSlice )

    if ReturnMesh:

        return xL, xU, nX, FileArray

    else:

        return
