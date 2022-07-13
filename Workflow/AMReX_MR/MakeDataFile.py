#!/usr/bin/env python3

import numpy as np
import os
from multiprocessing import Process, cpu_count

from UtilitiesModule import Overwrite, GetData, ChoosePlotFile, GetFileArray

def MakeDataFile( Field, PlotFileDirectory, DataFileDirectory, \
                  PlotFileBaseName, CoordinateSystem, \
                  SSi = -1, SSf = -1, nSS = -1, \
                  UsePhysicalUnits = True, \
                  MaxLevel = -1, Verbose = False ):

    """
    Generate a directory containing data files where each data file corresponds
    to a specific AMReX plotfile. Mesh data are stored in the header of each
    data file and can be accessed with the ReadHeader function.
    """

    if Verbose: print( '\nRunning MakeDataFile...\n' )

    if not DataFileDirectory[-1] == '/' : DataFileDirectory += '/'

    if Verbose:
        print( '\nDataFileDirectory: {:}\n'.format( DataFileDirectory ) )

    OW = Overwrite( DataFileDirectory )

    if OW:

        os.system( 'rm -rf {:}'.format( DataFileDirectory ) )
        os.system(  'mkdir {:}'.format( DataFileDirectory ) )

        if PlotFileDirectory[-1] != '/' : PlotFileDirectory += '/'

        if Verbose :
            print( '\nPlotFileDirectory: {:}\n'.format( PlotFileDirectory ) )

        PlotFileNameArray = GetFileArray( PlotFileDirectory, PlotFileBaseName )

        if SSi < 0: SSi = 0
        if SSf < 0: SSf = PlotFileNameArray.shape[0] - 1
        if nSS < 0: nSS = PlotFileNameArray.shape[0]

        PlotFileArray = []
        for i in range( nSS ):
            iSS = SSi + np.int64( ( SSf - SSi ) / ( nSS - 1 ) * i )
            PlotFile = str( PlotFileNameArray[iSS] )
            if PlotFile[-1] == '/' :
                PlotFileArray.append( PlotFile[0:-1] )
            else:
                PlotFileArray.append( PlotFile )
        PlotFileArray = np.array( PlotFileArray )

        TimeHeaderBase = '# Time []: '
        X1Base         = '# X1_C []: '
        X2Base         = '# X2_C []: '
        X3Base         = '# X3_C []: '
        dX1Base        = '# dX1  []: '
        dX2Base        = '# dX2  []: '
        dX3Base        = '# dX3  []: '
        if UsePhysicalUnits and CoordinateSystem == 'cartesian' :
            TimeHeaderBase = '# Time [ms]: '
            X1Base         = '# X1_C [km]: '
            X2Base         = '# X2_C [km]: '
            X3Base         = '# X3_C [km]: '
            dX1Base        = '# dX1  [km]: '
            dX2Base        = '# dX2  [km]: '
            dX3Base        = '# dX3  [km]: '
        elif UsePhysicalUnits and CoordinateSystem == 'spherical' :
            TimeHeaderBase = '# Time [ms]: '
            X1Base         = '# X1_C [km]: '
            X2Base         = '# X2_C [rad]: '
            X3Base         = '# X3_C [rad]: '
            dX1Base        = '# dX1  [km]: '
            dX2Base        = '# dX2  [rad]: '
            dX3Base        = '# dX3  [rad]: '

        Data, DataUnits, \
          X1, X2, X3, dX1, dX2, dX3, xL, xU, nX, Time \
            = GetData( PlotFileDirectory, PlotFileBaseName, Field, \
                       CoordinateSystem, UsePhysicalUnits, \
                       argv = [ 'a', PlotFileArray[0] ], \
                       MaxLevel = MaxLevel, \
                       ReturnTime = True, ReturnMesh = True )

        nDimsX = 1
        if( nX[1] > 1 ): nDimsX += 1
        if( nX[2] > 1 ): nDimsX += 1

        def loop( iLo, iHi ):

            for i in range( iLo, iHi ):

                PlotFile = PlotFileArray[i]

                DataFile = DataFileDirectory + PlotFile + '.dat'

                if Verbose:
                    print( 'Generating data file: {:} ({:}/{:})'.format \
                             ( DataFile, i+1, nSS ) )

                Data, DataUnits, \
                  X1, X2, X3, dX1, dX2, dX3, xL, xU, nX, Time \
                    = GetData( PlotFileDirectory, PlotFileBaseName, Field, \
                               CoordinateSystem, UsePhysicalUnits, \
                               argv = [ 'a', PlotFile ], \
                               MaxLevel = MaxLevel, \
                               ReturnTime = True, ReturnMesh = True )

                if   nDimsX == 1:
                    DataShape = '{:d}'.format( X1.shape[0] )
                elif nDimsX == 2:
                    DataShape = '{:d} {:d}'.format( X1.shape[0], X2.shape[0] )
                else:
                    exit( 'MakeDataFile not implemented for nDimsX > 2' )

                # Save multi-D array with np.savetxt. Taken from:
                # https://stackoverflow.com/questions/3685265/
                # how-to-write-a-multidimensional-array-to-a-text-file

                with open( DataFile, 'w' ) as FileOut:

                    FileOut.write( '# {:}\n'             .format( DataFile  ) )
                    FileOut.write( '# Array Shape: {:}\n'.format( DataShape ) )
                    FileOut.write( '# Data Units: {:}\n' .format( DataUnits ) )

                    TimeHeader = TimeHeaderBase + '{:.16e}\n'.format( Time )
                    FileOut.write( TimeHeader )

                    FileOut.write( X1Base )
                    for iX1 in range( X1.shape[0] ):
                        FileOut.write( str( X1[iX1] ) + ' ' )
                    FileOut.write( '\n' )

                    FileOut.write( X2Base )
                    for iX2 in range( X2.shape[0] ):
                        FileOut.write( str( X2[iX2] ) + ' ' )
                    FileOut.write( '\n' )

                    FileOut.write( X3Base )
                    for iX3 in range( X3.shape[0] ):
                        FileOut.write( str( X3[iX3] ) + ' ' )
                    FileOut.write( '\n' )

                    FileOut.write( dX1Base )
                    for iX1 in range( dX1.shape[0] ):
                        FileOut.write( str( dX1[iX1] ) + ' ' )
                    FileOut.write( '\n' )

                    FileOut.write( dX2Base )
                    for iX2 in range( dX2.shape[0] ):
                        FileOut.write( str( dX2[iX2] ) + ' ' )
                    FileOut.write( '\n' )

                    FileOut.write( dX3Base )
                    for iX3 in range( dX3.shape[0] ):
                        FileOut.write( str( dX3[iX3] ) + ' ' )
                    FileOut.write( '\n' )

                    np.savetxt( FileOut, Data )

                # end with open( DataFileName, 'w' ) as FileOut

            # end for i in range( nSS )

        # end of loop( iLo, iHi )

        # Adapted from:
        # https://www.benmather.info/post/2018-11-24-multiprocessing-in-python/

        nProc = cpu_count()

        processes = []

        for i in range( nProc ):
            iLo = np.int64( np.float64( i     ) / np.float64( nProc ) * nSS )
            iHi = np.int64( np.float64( i + 1 ) / np.float64( nProc ) * nSS )
            p = Process( target = loop, args = (iLo,iHi) )
            p.start()
            processes.append( p )

        [ p.join() for p in processes ]

    # end if OW

    return # end of MakeDataFile

def ReadHeader( DataFile ):

    f = open( DataFile )

    dum = f.readline()

    s = f.readline(); ind = s.find( ':' )+1
    DataShape = np.array( list( map( np.int64, s[ind:].split() ) ), np.int64 )

    s = f.readline(); ind = s.find( ':' )+1
    DataUnits = s[ind:]

    s = f.readline(); ind = s.find( ':' )+1
    Time = np.float64( s[ind:] )

    s = f.readline(); ind = s.find( ':' )+1
    X1_C = np.array( list( map( np.float64, s[ind:].split() ) ), np.float64 )

    s = f.readline(); ind = s.find( ':' )+1
    X2_C = np.array( list( map( np.float64, s[ind:].split() ) ), np.float64 )

    s = f.readline(); ind = s.find( ':' )+1
    X3_C = np.array( list( map( np.float64, s[ind:].split() ) ), np.float64 )

    s = f.readline(); ind = s.find( ':' )+1
    dX1 = np.array( list( map( np.float64, s[ind:].split() ) ), np.float64 )

    s = f.readline(); ind = s.find( ':' )+1
    dX2 = np.array( list( map( np.float64, s[ind:].split() ) ), np.float64 )

    s = f.readline(); ind = s.find( ':' )+1
    dX3 = np.array( list( map( np.float64, s[ind:].split() ) ), np.float64 )

    f.close()

    return DataShape, DataUnits, Time, X1_C, X2_C, X3_C, dX1, dX2, dX3
