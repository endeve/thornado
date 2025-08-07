import numpy as np
from sys import argv
from os import listdir, rename, mkdir, chdir
from shutil import copytree

# Based on Workflow/AMReX/Utilities/Files.py

DataDir = argv[1]
PlotBaseName = argv[2] + '.plt'
fileArray = []
PltfileArray = []

N = len( PlotBaseName )

chdir( DataDir )

fileArray \
  = listdir( './'  )

try:
    mkdir( 'ParaviewFiles' )
except:
    print('Sub-directory already exists.')

for iFile in range( len( fileArray ) ):
        
    sFile = fileArray[iFile]

    print( sFile )

    if ( sFile[0:N] == PlotBaseName and sFile[N+1].isdigit() ):
        
        copytree( './' + sFile, './ParaviewFiles/' + sFile )

        PltfileArray.append( sFile )

chdir( 'ParaviewFiles' )

for iFile  in range( len( PltfileArray ) ):

    pFile = PltfileArray[iFile]

    OldName = pFile
    
    print( str( pFile[N+1:] ) )

    NewName = 'plt' + str( pFile[N+1:] )

    rename( OldName, NewName )
