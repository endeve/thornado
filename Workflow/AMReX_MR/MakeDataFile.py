#!/usr/bin/env python3

import numpy as np
import os
from multiprocessing import Process, cpu_count

from UtilitiesModule import Overwrite, GetData, GetFileArray

def MakeDataFile( Field, PlotfileDirectory, DataDirectory, \
                  PlotfileBaseName, CoordinateSystem, \
                  SSi = -1, SSf = -1, nSS = -1, \
                  UsePhysicalUnits = True, \
                  MaxLevel = -1, \
                  forceChoiceD = False, owD = True, \
                  forceChoiceF = False, owF = True, \
                  Verbose = False ):

    """
    Generate a directory with the following structure for each plotfile
    (the example uses plotfile PlotfileBaseName.08675309)

    .DataDirectory/
    .DataDirectory/PlotfileNumbers.dat
    .DataDirectory/08675309/Time.dat
    .DataDirectory/08675309/X1.dat
    .DataDirectory/08675309/X2.dat
    .DataDirectory/08675309/X3.dat
    .DataDirectory/08675309/dX1.dat
    .DataDirectory/08675309/dX2.dat
    .DataDirectory/08675309/dX3.dat
    .DataDirectory/08675309/CF_D.dat
    .DataDirectory/08675309/PF_V1.dat
    .DataDirectory/08675309/<Your favorite field>.dat
    """

    print( '\nRunning MakeDataFile...' )

    if not DataDirectory[-1] == '/': DataDirectory += '/'

    print( '\nDataDirectory: {:}\n'.format( DataDirectory ) )

    ow = Overwrite( DataDirectory, ForceChoice = forceChoiceD, OW = owD )

    if ow:

        os.system( 'rm -rf {:}'.format( DataDirectory ) )
        os.system(  'mkdir {:}'.format( DataDirectory ) )

        PlotfileArray \
          = GetFileArray \
              ( PlotfileDirectory, PlotfileBaseName, \
                SSi = SSi, SSf = SSf, nSS = nSS )

        with open( DataDirectory + 'PlotfileNumbers.dat', 'w' ) as FileOut:
            for i in range( PlotfileArray.shape[0] ):
                FileOut.write( str(PlotfileArray[i][-8:]).zfill(8) + '\n' )

        if PlotfileDirectory[-1] != '/': PlotfileDirectory += '/'

        print( '\nPlotfileDirectory: {:}\n'.format( PlotfileDirectory ) )

        nSSS = nSS
        if nSS < 0: nSSS = PlotfileArray.shape[0]

        TimeUnits = '[]'
        X1Units   = '[]'
        X2Units   = '[]'
        X3Units   = '[]'
        if UsePhysicalUnits and CoordinateSystem == 'cartesian':
            TimeUnits = '[ms]'
            X1Units   = '[km]'
            X2Units   = '[km]'
            X3Units   = '[km]'
        elif UsePhysicalUnits and CoordinateSystem == 'spherical':
            TimeUnits = '[ms]'
            X1Units   = '[km]'
            X2Units   = '[rad]'
            X3Units   = '[rad]'

        TimeHeaderBase = '# Time {:}: '.format( TimeUnits )
        X1Base         = '# X1_C {:}: '.format( X1Units )
        X2Base         = '# X2_C {:}: '.format( X2Units )
        X3Base         = '# X3_C {:}: '.format( X3Units )
        dX1Base        = '# dX1  {:}: '.format( X1Units )
        dX2Base        = '# dX2  {:}: '.format( X2Units )
        dX3Base        = '# dX3  {:}: '.format( X3Units )

        printProcMem = False

        if printProcMem:
            import psutil
            process = psutil.Process( os.getpid() )
            print( 'mem: {:.3e} kB'.format \
                    ( process.memory_info().rss / 1024.0 ) )

        def loop( iLo, iHi ):

            for i in range( iLo, iHi ):

                if printProcMem:
                    print( 'mem: {:.3e} kB'.format \
                            ( process.memory_info().rss / 1024.0 ) )

                PlotfileNumber = PlotfileArray[i][-8:]

                FileDirectory = DataDirectory + PlotfileNumber + '/'

                os.system( 'mkdir {:}'.format( FileDirectory ) )

                TimeFile = FileDirectory + '{:}.dat'.format( 'Time' )
                X1File   = FileDirectory + '{:}.dat'.format( 'X1' )
                X2File   = FileDirectory + '{:}.dat'.format( 'X2' )
                X3File   = FileDirectory + '{:}.dat'.format( 'X3' )
                dX1File  = FileDirectory + '{:}.dat'.format( 'dX1' )
                dX2File  = FileDirectory + '{:}.dat'.format( 'dX2' )
                dX3File  = FileDirectory + '{:}.dat'.format( 'dX3' )
                DataFile = FileDirectory + '{:}.dat'.format( Field )

                if Verbose:
                    print( 'Generating data file: {:} ({:}/{:})'.format \
                             ( DataFile, i+1-iLo, iHi-iLo ) )

                Data, DataUnits, \
                  X1, X2, X3, dX1, dX2, dX3, xL, xH, nX, Time \
                    = GetData( PlotfileDirectory, PlotfileBaseName, Field, \
                               CoordinateSystem, UsePhysicalUnits, \
                               argv = [ 'a', PlotfileNumber ], \
                               MaxLevel = MaxLevel, \
                               SSi = SSi, SSf = SSf, nSS = nSS, \
                               ReturnTime = True, ReturnMesh = True )

                nDimsX = 1
                if( nX[1] > 1 ): nDimsX += 1
                if( nX[2] > 1 ): nDimsX += 1

                if   nDimsX == 1:
                    LoopShape = [ Data.shape[0], 1, 1 ]
                    DataShape = '{:d}' \
                                .format( Data.shape[0] )
                    Data = np.copy( Data[:,0  ,0  ] )
                    X1   = np.copy( X1  [:,0:1,0:1] )
                    X2   = np.copy( X2  [:,0:1,0:1] )
                    X3   = np.copy( X3  [:,0:1,0:1] )
                    dX1  = np.copy( dX1 [:,0:1,0:1] )
                    dX2  = np.copy( dX2 [:,0:1,0:1] )
                    dX3  = np.copy( dX3 [:,0:1,0:1] )
                elif nDimsX == 2:
                    LoopShape = [ Data.shape[0], Data.shape[1], 1 ]
                    DataShape = '{:d} {:d}' \
                                .format( Data.shape[0], Data.shape[1] )
                    Data = np.copy( Data[:,:,0  ] )
                    X1   = np.copy( X1  [:,:,0:1] )
                    X2   = np.copy( X2  [:,:,0:1] )
                    X3   = np.copy( X3  [:,:,0:1] )
                    dX1  = np.copy( dX1 [:,:,0:1] )
                    dX2  = np.copy( dX2 [:,:,0:1] )
                    dX3  = np.copy( dX3 [:,:,0:1] )
                else:
                    exit( 'MakeDataFile not implemented for nDimsX > 2' )


                # Save multi-D array with np.savetxt. Taken from:
                # https://stackoverflow.com/questions/3685265/
                # how-to-write-a-multidimensional-array-to-a-text-file

                with open( DataFile, 'w' ) as FileOut:

                    FileOut.write( '# {:}\n'             .format( DataFile  ) )
                    FileOut.write( '# Array Shape: {:}\n'.format( DataShape ) )
                    FileOut.write( '# Data Units: {:}\n' .format( DataUnits ) )
                    FileOut.write( '# Min. value: {:.16e}\n' \
                                   .format( Data.min() ) )
                    FileOut.write( '# Max. value: {:.16e}\n' \
                                   .format( Data.max() ) )

                    np.savetxt( FileOut, Data )

                # end with open( DataFile, 'w' ) as FileOut

                with open( TimeFile, 'w' ) as FileOut:

                    FileOut.write( '# {:}\n'              \
                                   .format( TimeFile  ) )
                    FileOut.write( '# Time Units: {:}\n'  \
                                   .format( TimeUnits ) )
                    FileOut.write( str( Time ) + '\n' )

                with open( X1File, 'w' ) as FileOut:

                    FileOut.write( '# {:}\n'.format( X1File  ) )
                    FileOut.write( '# X1_C {:}\n'.format( X1Units ) )

                    for iX1 in range( LoopShape[0] ):
                        for iX2 in range( LoopShape[1] ):
                            for iX3 in range( LoopShape[2] ):
                                FileOut.write( str( X1 [iX1,iX2,iX3] ) + ' ' )
                            FileOut.write( '\n' )
                        FileOut.write( '\n' )
                    FileOut.write( '\n' )

                with open( dX1File, 'w' ) as FileOut:

                    FileOut.write( '# {:}\n'.format( dX1File  ) )
                    FileOut.write( '# dX1 {:}\n'.format( X1Units ) )

                    for iX1 in range( LoopShape[0] ):
                        for iX2 in range( LoopShape[1] ):
                            for iX3 in range( LoopShape[2] ):
                                FileOut.write( str( dX1[iX1,iX2,iX3] ) + ' ' )
                            FileOut.write( '\n' )
                        FileOut.write( '\n' )
                    FileOut.write( '\n' )

                with open( X2File, 'w' ) as FileOut:

                    FileOut.write( '# {:}\n'.format( X2File  ) )
                    FileOut.write( '# X2_C {:}\n'.format( X2Units ) )

                    for iX1 in range( LoopShape[0] ):
                        for iX2 in range( LoopShape[1] ):
                            for iX3 in range( LoopShape[2] ):
                                FileOut.write( str( X2 [iX1,iX2,iX3] ) + ' ' )
                            FileOut.write( '\n' )
                        FileOut.write( '\n' )
                    FileOut.write( '\n' )

                with open( dX2File, 'w' ) as FileOut:

                    FileOut.write( '# {:}\n'.format( dX2File  ) )
                    FileOut.write( '# dX2 {:}\n'.format( X2Units ) )

                    for iX1 in range( LoopShape[0] ):
                        for iX2 in range( LoopShape[1] ):
                            for iX3 in range( LoopShape[2] ):
                                FileOut.write( str( dX2[iX1,iX2,iX3] ) + ' ' )
                            FileOut.write( '\n' )
                        FileOut.write( '\n' )
                    FileOut.write( '\n' )

                with open( X3File, 'w' ) as FileOut:

                    FileOut.write( '# {:}\n'.format( X3File  ) )
                    FileOut.write( '# X3_C {:}\n'.format( X3Units ) )

                    for iX1 in range( LoopShape[0] ):
                        for iX2 in range( LoopShape[1] ):
                            for iX3 in range( LoopShape[2] ):
                                FileOut.write( str( X3 [iX1,iX2,iX3] ) + ' ' )
                            FileOut.write( '\n' )
                        FileOut.write( '\n' )
                    FileOut.write( '\n' )

                with open( dX3File, 'w' ) as FileOut:

                    FileOut.write( '# {:}\n'.format( dX3File  ) )
                    FileOut.write( '# dX3 {:}\n'.format( X3Units ) )

                    for iX1 in range( LoopShape[0] ):
                        for iX2 in range( LoopShape[1] ):
                            for iX3 in range( LoopShape[2] ):
                                FileOut.write( str( dX3[iX1,iX2,iX3] ) + ' ' )
                            FileOut.write( '\n' )
                        FileOut.write( '\n' )
                    FileOut.write( '\n' )

            # end for i in range( iLo, iHi )

        # end of loop( iLo, iHi )

        # Adapted from:
        # https://www.benmather.info/post/2018-11-24-multiprocessing-in-python/

        nProc = max( 1, cpu_count() // 2 )

        print( 'Generating {:} with {:} processes...\n'.format \
             ( DataDirectory, nProc ) )

        if nProc > 1:

          processes = []

          for i in range( nProc ):
              iLo = np.int64( np.float64( i     ) / np.float64( nProc ) * nSSS )
              iHi = np.int64( np.float64( i + 1 ) / np.float64( nProc ) * nSSS )
              p = Process( target = loop, args = (iLo,iHi) )
              p.start()
              processes.append( p )

          [ p.join() for p in processes ]

    else:

        PlotfileArray \
          = np.loadtxt( DataDirectory + 'PlotfileNumbers.dat', dtype = str )

        File = DataDirectory + PlotfileArray[0] + '/{:}.dat'.format( Field )

        oww = Overwrite( File, ForceChoice = forceChoiceF, OW = owF )

        if oww:

            print( '\nPlotfileDirectory: {:}\n'.format( PlotfileDirectory ) )

            nSSS = nSS
            if nSS < 0: nSSS = PlotfileArray.shape[0]

            printProcMem = False

            if printProcMem:
                import psutil
                process = psutil.Process( os.getpid() )
                print( 'mem: {:.3e} kB'.format \
                        ( process.memory_info().rss / 1024.0 ) )

            def loop( iLo, iHi ):

                for i in range( iLo, iHi ):

                    if printProcMem:
                        print( 'mem: {:.3e} kB'.format \
                                ( process.memory_info().rss / 1024.0 ) )

                    PlotfileNumber = PlotfileArray[i]

                    FileDirectory = DataDirectory + PlotfileNumber + '/'

                    DataFile = FileDirectory + '{:}.dat'.format( Field )

                    if Verbose:
                        print( 'Generating data file: {:} ({:}/{:})'.format \
                                 ( DataFile, i+1-iLo, iHi-iLo ) )

                    Data, DataUnits, \
                      X1, X2, X3, dX1, dX2, dX3, xL, xH, nX \
                        = GetData( PlotfileDirectory, PlotfileBaseName, Field, \
                                   CoordinateSystem, UsePhysicalUnits, \
                                   argv = [ 'a', PlotfileNumber ], \
                                   MaxLevel = MaxLevel, \
                                   SSi = SSi, SSf = SSf, nSS = nSS, \
                                   ReturnTime = False, ReturnMesh = True )

                    nDimsX = 1
                    if( nX[1] > 1 ): nDimsX += 1
                    if( nX[2] > 1 ): nDimsX += 1

                    if   nDimsX == 1:
                        LoopShape = [ Data.shape[0], 1, 1 ]
                        DataShape = '{:d}' \
                                    .format( Data.shape[0] )
                        Data = np.copy( Data[:,0  ,0  ] )
                    elif nDimsX == 2:
                        LoopShape = [ Data.shape[0], Data.shape[1], 1 ]
                        DataShape = '{:d} {:d}' \
                                    .format( Data.shape[0], Data.shape[1] )
                        Data = np.copy( Data[:,:,0  ] )
                    else:
                        exit( 'MakeDataFile not implemented for nDimsX > 2' )


                    # Save multi-D array with np.savetxt. Taken from:
                    # https://stackoverflow.com/questions/3685265/
                    # how-to-write-a-multidimensional-array-to-a-text-file

                    with open( DataFile, 'w' ) as FileOut:

                        FileOut.write( '# {:}\n'             .format( DataFile  ) )
                        FileOut.write( '# Array Shape: {:}\n'.format( DataShape ) )
                        FileOut.write( '# Data Units: {:}\n' .format( DataUnits ) )
                        FileOut.write( '# Min. value: {:.16e}\n' \
                                       .format( Data.min() ) )
                        FileOut.write( '# Max. value: {:.16e}\n' \
                                       .format( Data.max() ) )

                        np.savetxt( FileOut, Data )

                    # end with open( DataFile, 'w' ) as FileOut

                # end for i in range( iLo, iHi )

            # end of loop( iLo, iHi )

            # Adapted from:
            # https://www.benmather.info/post/2018-11-24-multiprocessing-in-python/

            nProc = max( 1, cpu_count() // 2 )

            print( 'Generating {:} with {:} processes...\n'.format \
                 ( DataDirectory, nProc ) )

            if nProc > 1:

              processes = []

              for i in range( nProc ):
                  iLo = np.int64( np.float64( i     ) / np.float64( nProc ) * nSSS )
                  iHi = np.int64( np.float64( i + 1 ) / np.float64( nProc ) * nSSS )
                  p = Process( target = loop, args = (iLo,iHi) )
                  p.start()
                  processes.append( p )

              [ p.join() for p in processes ]

            else:

              loop( 0, nSSS )
        # END if oww

    # END if ow

    PlotfileArray \
      = np.loadtxt( DataDirectory + 'PlotfileNumbers.dat', dtype = str )

    return PlotfileArray
# END of MakeDataFile

def ReadHeader( DataFile ):

    f = open( DataFile )

    dum = f.readline()

    s = f.readline(); ind = s.find( ':' )+1
    DataShape = np.array( list( map( np.int64, s[ind:].split() ) ), np.int64 )

    s = f.readline(); ind = s.find( ':' )+1
    DataUnits = s[ind:]

    s = f.readline(); ind = s.find( ':' )+1
    MinVal = np.float64( s[ind:] )

    s = f.readline(); ind = s.find( ':' )+1
    MaxVal = np.float64( s[ind:] )

    f.close()

    return DataShape, DataUnits, MinVal, MaxVal
# END of ReadHeader
