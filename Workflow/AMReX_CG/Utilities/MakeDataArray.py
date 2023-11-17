#!/usr/bin/env python3

import numpy as np
import os
from os.path import isdir, isfile
from multiprocessing import Process, cpu_count

import GlobalVariables.Settings as gvS
import GlobalVariables.Units    as gvU

from Utilities.Files        import Overwrite, CheckForField, CleanField
from Utilities.GetFrameData import GetFrameData
from Utilities.Mesh         import CreateLocations


#=============================================#
#   Included Routines
#
#   MakeProbelmDataDirectory
#   MakeDataFrameDirectories
#   MakeDataFrameDirectoriesLoop
#   CreateFrameDirectory
#
#
#
#  Code Structure
#
#   MakeProbelmDataDirectory             ( Queries for overwrite states.   )
#   \=> MakeDataFrameDirectories         ( Divides work amoung available processes.)
#       \=> MakeDataFrameDirectoriesLoop ( Applies overwriteF.  )
#           \=> CreateFrameDirectory     ( Creates and fills individual frame directory.)
#
#
#
#  Data Directory Structure
#
#   .DataDirectory/                          ( Usually named for the problem. )
#   .DataDirectory/FrameNumber/              ( FrameNumber is an 8 digit integer. )
#   .DataDirectory/FrameNumber/X1.dat        ( Mesh information. )
#   .DataDirectory/FrameNumber/X2.dat        ( Mesh information. )
#   .DataDirectory/FrameNumber/X3.dat        ( Mesh information. )
#   .DataDirectory/FrameNumber/dX1.dat       ( Mesh information. )
#   .DataDirectory/FrameNumber/dX2.dat       ( Mesh information. )
#   .DataDirectory/FrameNumber/dX3.dat       ( Mesh information. )
#   .DataDirectory/FrameNumber/FieldName.dat ( Field data.       )
#
#==============================================================#


 #=============================================#
#                                               #
#   MakeProbelmDataDirectory                    #
#                                               #
 #=============================================#
def MakeProbelmDataDirectory( FileNumberArray,          \
                              PlotDirectory,            \
                              PlotBaseName,             \
                              Field,                    \
                              DataDirectory,            \
                              forceChoiceD     = False, \
                              overwriteD       = True,  \
                              forceChoiceF     = False, \
                              overwriteF       = True,  \
                              SaveTime         = True,  \
                              SaveX            = True,  \
                              SavedX           = True   ):
                     
    """
    Generate a directory with the following structure for each plotfile
    (the example uses plotfile PlotBaseName.08675309)

    .DataDirectory/
    .DataDirectory/PlotfileNumbers.dat
    .DataDirectory/08675309/Time.dat
    .DataDirectory/08675309/X1.dat
    .DataDirectory/08675309/X2.dat
    .DataDirectory/08675309/X3.dat
    .DataDirectory/08675309/CF_D.dat
    .DataDirectory/08675309/PF_V1.dat
    .DataDirectory/08675309/<Your favorite field>.dat
    """

    print( '\n  Running MakeDataFile' )
    
    print(   '  --------------------' )

    if DataDirectory[-1] != '/': DataDirectory += '/'
    if PlotDirectory[-1] != '/': PlotDirectory += '/'

    print( '\n  DataDirectory: {:}\n'.format( DataDirectory ) )
    print( '\n  PlotDirectory: {:}\n'.format( PlotDirectory ) )
    
    
    overwriteD = Overwrite( DataDirectory, ForceChoice = forceChoiceD, OW = overwriteD )

    

    nProc = 1 #max( 1, cpu_count() // 2 )

    if overwriteD:
    
        
        os.system( 'rm -rf {:}'.format( DataDirectory ) )
        os.system(  'mkdir {:}'.format( DataDirectory ) )
    
                    
        print( '  Generating {:} with {:} processes...\n'.format \
             ( DataDirectory, nProc ) )

            

    else: # overwriteD == False


        #   Check to see if any of the existing data directories
        # contain data related to the specified field.
        PathDataDirectory, Flag = CheckForField(Field,          \
                                                DataDirectory,  \
                                                FileNumberArray     )
        
        
        #   If data related to specified field exists, ask user
        # if they wish it to be overwritten.
        if Flag:
            overwriteF = Overwrite( PathDataDirectory,         \
                                    ForceChoice = forceChoiceF,\
                                    OW = overwriteF            )

        
        if overwriteF:
        
            #   If user wishes to overwrite existing data,
            # this cleans the existing data directories.
            CleanField( Field,          \
                        DataDirectory,  \
                        FileNumberArray     )
        
    
    
    #   Creates a data directory with sub-directories containing
    # data needed to produce plots.
    MakeDataFrameDirectories(   FileNumberArray,    \
                                PlotDirectory,      \
                                PlotBaseName,       \
                                Field,              \
                                DataDirectory,      \
                                SaveTime = SaveTime,\
                                SaveX = SaveX,      \
                                SavedX = SavedX,    \
                                nProcs = nProc,     \
                                owF = overwriteF,   )

    return





 #=============================================#
#                                               #
#   MakeDataFrameDirectories                    #
#                                               #
 #=============================================#
def MakeDataFrameDirectories( FileNumberArray,          \
                              PlotDirectory,            \
                              PlotBaseName,             \
                              Field,                    \
                              DataDirectory,            \
                              nProcs           = 1,     \
                              owF              = True,  \
                              SaveTime         = True,  \
                              SaveX            = True,  \
                              SavedX           = True   ):

    nFiles = FileNumberArray.shape[0]

    if nProcs > 1:
        processes = []
        for i in range( nProcs ):

            iLo = np.int64( np.float64( i     ) \
                  / np.float64( nProcs ) * nFiles )
            iHi = np.int64( np.float64( i + 1 ) \
                  / np.float64( nProcs ) * nFiles )
            
        
            p = Process( target=MakeDataFrameDirectoriesLoop,  \
                         args = (FileNumberArray[iLo:iHi],   \
                                 PlotDirectory,     \
                                 PlotBaseName,      \
                                 Field,             \
                                 DataDirectory,     \
                                 SaveTime,          \
                                 SaveX,             \
                                 SavedX,            \
                                 owF,               \
                                 gvS.Verbose,       )          )
            
            
            p.start()
            processes.append( p )
        
        [ p.join() for p in processes ]


              # Ensure all files were created

        fileArray = []
        for i in range( nFiles ):

            PlotfileNumber = str(FileNumberArray[i])
            FileDirectory = DataDirectory + PlotfileNumber + '/'

            TimeFile = FileDirectory + '{:}.dat'.format( 'Time' )
            X1File   = FileDirectory + '{:}.dat'.format( 'X1' )
            X2File   = FileDirectory + '{:}.dat'.format( 'X2' )
            X3File   = FileDirectory + '{:}.dat'.format( 'X3' )
            dX1File  = FileDirectory + '{:}.dat'.format( 'dX1' )
            dX2File  = FileDirectory + '{:}.dat'.format( 'dX2' )
            dX3File  = FileDirectory + '{:}.dat'.format( 'dX3' )
            DataFile = FileDirectory + '{:}.dat'.format( Field )

            if not isfile( TimeFile ) or not isfile( X1File ) \
                or not isfile( X2File ) or not isfile( X3File ) \
                or not isfile( dX1File ) or not isfile( dX2File ) \
                or not isfile( dX3File ) or not isfile( DataFile ):
                    fileArray.append( i )

        if len( fileArray ) != 0:
            MakeDataFrameDirectoriesLoop( fileArray,      \
                                          PlotDirectory,  \
                                          PlotBaseName,   \
                                          Field,          \
                                          DataDirectory,  \
                                          SaveTime,       \
                                          SaveX,          \
                                          SavedX,         \
                                          owF,            \
                                          gvS.Verbose,    )


    else:
        MakeDataFrameDirectoriesLoop( FileNumberArray,    \
                                      PlotDirectory,      \
                                      PlotBaseName,       \
                                      Field,              \
                                      DataDirectory,      \
                                      SaveTime,           \
                                      SaveX,              \
                                      SavedX,             \
                                      owF,                \
                                      gvS.Verbose,        )










 #=============================================#
#                                               #
#   MakeDataFrameDirectoriesLoop                #
#                                               #
 #=============================================#
def MakeDataFrameDirectoriesLoop( FileNumberArray,  \
                                  PlotDirectory,    \
                                  PlotBaseName,     \
                                  Field,            \
                                  DataDirectory,    \
                                  SaveTime,         \
                                  SaveX,            \
                                  SavedX,           \
                                  owF,              \
                                  Verbose           ):

    

#    NumPltFiles = FileNumberArray.shape[0]
    NumPltFiles = len(FileNumberArray)
    
    for i in range(NumPltFiles):

    
        printProcMem = False
        if printProcMem:
            print( 'mem: {:.3e} kB'.format \
                    ( process.memory_info().rss / 1024.0 ) )

        PlotFileNumber = FileNumberArray[i]
        PathDataDirectory = DataDirectory + str(PlotFileNumber) + '/'
        PathPlotDirectory = PlotDirectory       \
                          + PlotBaseName    \
                          + '{:}'.format( str(PlotFileNumber).zfill(8) )
        
        if not isdir(PathDataDirectory):
        
            if Verbose:
                print( 'Generating data directory: {:} ({:}/{:})'.format \
                    ( PathDataDirectory, i+1, NumPltFiles ) )
                    
            CreateFrameDirectory( PathDataDirectory, \
                                  PathPlotDirectory, \
                                  PlotFileNumber,    \
                                  Field,             \
                                  SaveTime,          \
                                  SaveX,             \
                                  SavedX             )
                             
        elif isdir(PathDataDirectory) and owF:
        
            PathFieldDirectory = PathDataDirectory         \
                               + '{:}.dat'.format( Field )

            os.system( 'rm -rf {:}'.format( PathFieldDirectory ) )
            
            
            
            if Verbose:
                print( 'Generating data directory: {:} ({:}/{:})'.format \
                    ( PathDataDirectory, i+1, NumPltFiles ) )
                    
            CreateFrameDirectory( PathDataDirectory, \
                                  PathPlotDirectory, \
                                  PlotFileNumber,    \
                                  Field,             \
                                  SaveTime,          \
                                  SaveX,             \
                                  SavedX             )
        
    return









 #=============================================#
#                                               #
#   CreateFrameDirectory                        #
#                                               #
 #=============================================#
def CreateFrameDirectory( PathDataDirectory,             \
                          PathPlotDirectory,             \
                          PlotFileNumber,                \
                          Field,                         \
                          SaveTime         = True,       \
                          SaveX            = True,       \
                          SavedX           = True        ):

    if not isdir(PathDataDirectory):
        os.system( 'mkdir {:}'.format(PathDataDirectory) )

    DataFile = PathDataDirectory + '{:}.dat'.format( Field )

    if SaveTime:
        TimeFile = PathDataDirectory + '{:}.dat'.format( 'Time' )
    if SaveX:
        X1File   = PathDataDirectory + '{:}.dat'.format( 'X1' )
        X2File   = PathDataDirectory + '{:}.dat'.format( 'X2' )
        X3File   = PathDataDirectory + '{:}.dat'.format( 'X3' )
    if SavedX:
        dX1File  = PathDataDirectory + '{:}.dat'.format( 'dX1' )
        dX2File  = PathDataDirectory + '{:}.dat'.format( 'dX2' )
        dX3File  = PathDataDirectory + '{:}.dat'.format( 'dX3' )

    

    X1, X2, X3, dX1, dX2, dX3, xL, xH               \
        = CreateLocations(  PathPlotDirectory,      \
                            TypeIn = "Leaf"      )

    Data, DataUnits, Time \
        = GetFrameData( PathPlotDirectory,   \
                        Field,               \
                        X1, X2, X3,          \
                        dX1, dX2, dX3,       \
                        SaveTime             )


    nX1 = X1.shape[0]
    nX2 = X2.shape[0]
    nX3 = X3.shape[0]


    nDims = 1
    if nX2 > 1: nDims += 1
    if nX3 > 1: nDims += 1

    if   nDims == 1:
        LoopShape = [ Data.shape[0], 1, 1 ]
        DataShape = '{:d}' \
                    .format( Data.shape[0] )

    else:
        exit( 'MakeDataFile not implemented for nDims > 1' )



    if not isfile( DataFile ):

        # Save multi-D array with np.savetxt. Taken from:
        # https://stackoverflow.com/questions/3685265/
        # how-to-write-a-multidimensional-array-to-a-text-file

        with open( DataFile, 'w' ) as FileOut:

            FileOut.write( '# {:}\n'              \
                           .format( DataFile  ) )
            FileOut.write( '# Array Shape: {:}\n' \
                           .format( DataShape ) )
            FileOut.write( '# Data Units: {:}\n'  \
                           .format( DataUnits ) )
            FileOut.write( '# Min. value: {:.16e}\n' \
                           .format( Data.min() ) )
            FileOut.write( '# Max. value: {:.16e}\n' \
                           .format( Data.max() ) )

            np.savetxt( FileOut, Data )

        # end with open( DataFile, 'w' ) as FileOut

    # end if not isfile( DataFile )

    if SaveTime:
        if not isfile( TimeFile ):

            with open( TimeFile, 'w' ) as FileOut:

                FileOut.write( '# {:}\n'              \
                               .format( TimeFile  ) )
                FileOut.write( '# Time Units: {:}\n'  \
                               .format( gvU.TimeUnits ) )
                FileOut.write( str( Time ) + '\n' )




    if SaveX:
        if not isfile( X1File ):

            with open( X1File, 'w' ) as FileOut:

                FileOut.write( '# {:}\n'.format( X1File  ) )
                FileOut.write( '# X1_C [{:}]\n'.format( gvU.X1Units ) )

                for iX1 in range( LoopShape[0] ):
                    FileOut.write( str( X1 [iX1] ) + ' ' )
                FileOut.write( '\n' )


        if not isfile( X2File ):

            with open( X2File, 'w' ) as FileOut:

                FileOut.write( '# {:}\n'.format( X2File  ) )
                FileOut.write( '# X2_C [{:}]\n'.format( gvU.X2Units ) )

                for iX2 in range( LoopShape[1] ):
                    FileOut.write( str( X2 [iX2] ) + ' ' )
                FileOut.write( '\n' )


        if not isfile( X3File ):

            with open( X3File, 'w' ) as FileOut:

                FileOut.write( '# {:}\n'.format( X3File  ) )
                FileOut.write( '# X3_C [{:}]\n'.format( gvU.X3Units ) )

                for iX3 in range( LoopShape[2] ):
                    FileOut.write \
                      ( str( X3 [iX3] ) + ' ' )
                FileOut.write( '\n' )




    if SavedX:
        if not isfile( dX1File ):

            with open( dX1File, 'w' ) as FileOut:

                FileOut.write( '# {:}\n'.format( dX1File  ) )
                FileOut.write( '# dX1 [{:}]\n'.format( gvU.X1Units ) )

                for iX1 in range( LoopShape[0] ):
                    FileOut.write( str( dX1 [iX1] ) + ' ' )
                FileOut.write( '\n' )


        if not isfile( dX2File ):

            with open( dX2File, 'w' ) as FileOut:

                FileOut.write( '# {:}\n'.format( dX2File  ) )
                FileOut.write( '# dX2 [{:}]\n'.format( gvU.X2Units ) )

                for iX2 in range( LoopShape[1] ):
                    FileOut.write( str( dX2 [iX2] ) + ' ' )
                FileOut.write( '\n' )


        if not isfile( dX3File ):

            with open( dX3File, 'w' ) as FileOut:

                FileOut.write( '# {:}\n'.format( dX3File  ) )
                FileOut.write( '# dX3 [{:}]\n'.format( gvU.X3Units ) )

                for iX3 in range( LoopShape[2] ):
                    FileOut.write( str( dX3 [iX3] ) + ' ' )
                FileOut.write( '\n' )
