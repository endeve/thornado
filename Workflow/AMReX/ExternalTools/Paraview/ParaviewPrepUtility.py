from os import chdir, listdir, mkdir
from os import symlink, rename, system
from sys import argv
from shutil import copytree

# Based on Workflow/AMReX/Utilities/Files.py

ScratchDir = argv[1] # The path to the directory containing multiple datasets.
DataDir = argv[2] # The directory containing a single dataset.
PlotBaseName = argv[3] + '.plt'
fileArray = []
PltfileArray = []

N = len( PlotBaseName )

chdir( ScratchDir + '/' + DataDir )

fileArray \
  = listdir( './'  )

# Make directories for data symlinks and metadata.

try:
    mkdir( 'ParaviewFiles' )
except:
    print( 'Paraview ub-directory already exists.' )

try:
    mkdir( 'ParaviewFiles/Metadata' )
except:
    print( 'Metadata sub-directory already exists.' )

# Make symlinks and put in ParaviewFiles directory.

for iFile in range( len( fileArray ) ):
        
    oFile = fileArray[iFile]

    print( str( oFile[N+1:] ) )

    NewName = 'plt' + str( oFile[N+1:] )

    if ( oFile[0:N] == PlotBaseName and oFile[N+1].isdigit() ):
       
        try:
            symlink( './' + oFile, './ParaviewFiles/' + NewName, target_is_directory = True )
            PltfileArray.append( oFile )
        except:
                print('Symlink already exists.')

# tar up full directory

chdir( ScratchDir )
try:
    system( 'tar -czvf ' + DataDir + '_ParaviewPlt.tar.gz' + ' ' + DataDir + '/ParaviewFiles' + ' ' + DataDir + '/*.plt*' )
except:
    print( 'tar has already been performed.' )
