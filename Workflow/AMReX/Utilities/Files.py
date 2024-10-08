#!/usr/bin/env python3

import numpy as np
import GlobalVariables.Settings as gvS
from charade import detect

#=============================================#
#   Included Routines
#
#   GetFileNumberArray
#   convert
#   ChoosePlotFile
#   Overwrite
#   CheckForField
#   CleanField
#
#=============================================#


#=============================================#
#                                             #
#   GetFileNumberArray                        #
#                                             #
#=============================================#
def GetFileNumberArray \
      ( PlotDirectory, PlotBaseName, \
        SSi               = -1, \
        SSf               = -1, \
        PlotEvery         = 1 , \
        snapshotIndexList = [] ):

    from os import listdir

    N = len( PlotBaseName )

    # Create list of all files in plotDirectory
    fileArray \
      = np.sort( np.array( \
          [ convert( file ) for file in listdir( PlotDirectory ) ] ) )

    # Filter file list to only those that start with plotBaseName
    fileList = []
    for iFile in range( fileArray.shape[0] ):

        sFile = fileArray[iFile]

        if ( sFile[0:N] == PlotBaseName and sFile[N+1].isdigit() ) :
            fileList.append( np.int64( sFile[N:N+8] ) )
        elif ( sFile[-3:] == '.h5' ) :
            fileList.append( int(sFile[-9:-3]) )

    

    fileNumbers = np.array( fileList )
    numNumbers = fileNumbers.shape[0]

    # Warn if list is empty.
    if not fileNumbers.shape[0] > 0:

        msg = '\n>>>No files found.\n'
        msg += '>>>Double check the path: {:}\n'.format( PlotDirectory )
        msg += '>>>Double check the PlotBaseName: {:}\n' \
               .format( PlotBaseName )
        msg += '>>> Is it plt_ or just plt?\n'

        assert ( fileNumbers.shape[0] > 0 ), msg

    # Sort Numbers
    fileNumbers.sort()

    # Filter File List
    if SSi < 0: SSi = np.int64( fileNumbers[0] )
    if SSf < 0: SSf = np.int64( fileNumbers[-1] )


    if SSf < SSi:
        msg = '\n>>> Final frame comes before initial frame. \n'
        msg += '>>> Check SSf > SSi.'
        assert ( SSf > SSi ), msg

    for SSi_index in range( numNumbers ):
        if ( fileNumbers[SSi_index] == SSi ): break

    for SSf_index in range( numNumbers-1, SSi_index-2,-1):
        if ( fileNumbers[SSf_index] == SSf ): break

    if not SSi_index <= numNumbers-1:
        msg = '\n>>> SSi does not correspond to a plot data directory.\n'
        assert( SSi_index < numNumbers-1 ), msg

    if not SSf_index >= SSi_index:
        msg = '\n>>> SSf does not correspond to a plot data directory.\n'
        assert( SSf_index >= SSi_index ), msg

    nSS = SSf_index - SSi_index + 1


    #   Filter file list to specified range
    fileNumbersLimited = []
    for i in range( nSS ):
        fileNumbersLimited.append(fileNumbers[SSi_index+i])


    if ( len( snapshotIndexList ) == 0 ) :
        # Filter file number list by PlotEvery value.
        Here = 0
        There = len(fileNumbersLimited)
        fileNumbersFiltered \
          = np.array( fileNumbersLimited[Here:There:PlotEvery] )
    else :
        # Filter file number list by user defined snapshotIndexList
        fileNumbersFiltered \
          = np.empty( len( snapshotIndexList ), dtype = np.int64 )
        for i in range( len( snapshotIndexList ) ):
            fileNumbersFiltered[i] = fileNumbers[snapshotIndexList[i]]


    return fileNumbersFiltered
# END GetFileNumberArray


#=============================================#
#                                             #
#   convert                                   #
#                                             #
#=============================================#
def convert( s ):

    """from https://www.geeksforgeeks.org/python-character-encoding/"""

    # if in the charade instance
    if isinstance(s, str):
        s = s.encode()

    # retrieving the encoding information
    # from the detect() output
    encoding = detect(s)['encoding']

    if encoding == 'utf-8':
        return s.decode()
    else:
        return s.decode(encoding)


#=============================================#
#                                             #
#   ChoosePlotFile                            #
#                                             #
#=============================================#
def ChoosePlotFile( FileNumberArray,            \
                    PlotBaseName = 'plt',       \
                    argv = [ 'a' ],             \
                    DataType = "AMReX",         \
                    Verbose = False ):


    if len( argv ) == 1:

        # Get last plotfile in directory
        File = PlotBaseName + '{:}'.format( str(FileNumberArray[-1]).zfill(8) )

    elif len( argv ) == 2:
        if argv[1][0].isalpha():

            File = argv[1]

        else:
            if DataType.lower() == 'amrex':
                File = PlotBaseName + '{:}'.format( argv[1].zfill(8) )
            else:
                File = PlotBaseName[:-4]+'_{:}'.format( argv[1].zfill(6) )
                
    else:

        n = len( argv )

        msg = 'len( argv ) must be > 0 and < 3: len( argv ) = {:d}'.format( n )

        arg = ( n > 0 ) & ( n < 3 )
        assert arg, msg

    # Remove "/" at end of filename, if present
    if File[-1] == '/' : File = np.copy( File[:-1] )

    if Verbose: print( File )

    return File


#=============================================#
#                                             #
#   Overwrite                                 #
#                                             #
#=============================================#
def Overwrite( FileOrDirName, ForceChoice = False, OW = False ):

    if ForceChoice: return OW

    from os.path import isfile, isdir

    OW = True

    if ( isfile( FileOrDirName ) or isdir( FileOrDirName ) ):

        if ( isdir( FileOrDirName ) and FileOrDirName[-1] != '/' ):
            FileOrDirName += '/'

        YN = input( '  {:} exists. Overwrite? (y/N): '.format( FileOrDirName ) )

        if YN == 'Y' or YN == 'y' :
            print( '  Overwriting' )
            OW = True
        else:
            print( '  Not overwriting' )
            OW = False

    return OW


#=============================================#
#                                             #
#   CheckForField                             #
#                                             #
#=============================================#
def CheckForField(  Field,          \
                    DataDirectory,  \
                    FileNumberArray ):

    from os.path import isfile

    if not DataDirectory[-1] == '/': DataDirectory += '/'

    Found_Flag = False
    PathToFirstField = ''
    for i in range(FileNumberArray.shape[0]):

        FilePath = DataDirectory + str(FileNumberArray[i]) + '/' + '{:}.dat'.format(Field)

        if isfile( FilePath ):
            Found_Flag = True
            PathToFirstField = FilePath
            break

    return PathToFirstField, Found_Flag


#=============================================#
#                                             #
#   CleanField                                #
#                                             #
#=============================================#
def CleanField( Field,          \
                DataDirectory,  \
                FileNumberArray ):

    import os
    from os      import scandir
    from os.path import isfile

    # Create list of all directories in DataDirectory
    DirectoryArray = [f.path for f in scandir(DataDirectory) if f.is_dir()]


    for i in range(len(DirectoryArray)):
        FilePath = DirectoryArray[i] + '/' + '{:}.dat'.format(Field)

        if isfile(FilePath):
            os.system('rm {:}'.format(FilePath))

    return





#=============================================#
#                                             #
#   CleanField                                #
#                                             #
#=============================================#
def CreatePlotPath( Directory,      \
                    BaseName,       \
                    FileNumber,     \
                    DataType        ):

    if DataType.lower() == 'amrex':
        Path = Directory        \
             + BaseName         \
             + '{:}'.format( str(FileNumber).zfill(8) )
    elif DataType.lower() == 'native':
        Path = Directory        \
             + BaseName[:-4]    \
             + '_{:}'.format( str(FileNumber).zfill(6) )
    else:
        print("Error in CreatePlotPath. Bad DataType. Input: ",DataType)

    return Path
