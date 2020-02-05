#!/usr/local/bin/python3

"""
Given a directory with 2D AMReX plot-files
in spherical-polar coordinates, create a data
file containing all the values of a specified
variable for each time-step.
"""

import yt
import numpy as np
import subprocess
from os import listdir
from os.path import isfile
from sys import argv, exit

yt.funcs.mylog.setLevel(0) # Suppress initial yt output to screen

# --- Get user's HOME directory ---
HOME = subprocess.check_output( ["echo $HOME"], shell = True)
HOME = HOME[:-1].decode( "utf-8" ) + '/'

############# User input #############
DataDirectory = HOME + '/'

TimeFileName     = 'MovieTime.dat'
UsePhysicalUnits = True # Are you using physical units?
Relativistic     = True # Are you plotting results from relativistic hydro?

PlotFileBaseName = 'thornado'

WriteExtras = False # Set to true if generating data on
                    # external machine (Summit, ACF, etc.)

############# End of user input #############

if( UsePhysicalUnits ):
    c = 2.99792458e10
    Centimeter = 1.0e5 # Centimeters per kilometer
else:
    c = 1.0
    Centimeter = 1.0

# Get array of all plot-files
FileArray = np.sort(np.array( [ file for file in listdir( DataDirectory ) ] ) )
FileList = []
for iFile in range( FileArray.shape[0] ):
    sFile = FileArray[iFile]
    if( sFile[0:len(PlotFileBaseName)+1] == PlotFileBaseName + '_' \
          and sFile[len(PlotFileBaseName)+1].isdigit() ):
        FileList.append( sFile )
FileArray = np.array( FileList )

# Get some general info about the computational domain
ds = yt.load( '{:}'.format( DataDirectory + FileArray[0] ) )
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

def MakeDataFile( Field, DataFileName ):

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

            if( Field == 'Entropy' ):

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
                  = CoveringGrid['PF_V1'   ].to_ndarray()[:,:,0] * Centimeter
                PF_V2 \
                  = CoveringGrid['PF_V2'   ].to_ndarray()[:,:,0] * Centimeter
                PF_V3 \
                  = CoveringGrid['PF_V3'   ].to_ndarray()[:,:,0] * Centimeter
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
                if( Relativistic ): # Plot h/c^2
                    Data[i] = ( c**2 + ( PF_E + AF_P ) / PF_D ) / c**2

            else:

                Data[i] \
                  = CoveringGrid[Field].to_ndarray()[:,:,0]

                if( Field[-2] == 'V' ):
                    Data[i] *= Centimeter

            Time[i] = ds.current_time

        # Save multi-D array with np.savetxt. Taken from:
        # https://stackoverflow.com/questions/3685265/how-to-write-a-multidimensional-array-to-a-text-file

        with open( DataFileName, 'w' ) as FileOut:
            FileOut.write( '# Array shape: {0}\n'.format( Data.shape ) )

            # Iterating through an n-dimensional array produces slices along
            # the last axis. This is equivalent to Data[i] in this case
            for TimeSlice in Data:
                np.savetxt( FileOut, TimeSlice )
                FileOut.write( '# New slice\n' )

        np.savetxt( TimeFileName, Time )

    return
