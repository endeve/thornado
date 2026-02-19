from os import chdir, listdir, mkdir, rename
from sys import argv

# Renames old AMReX plot and checkpoint files from
# ProblemName.plt# and ProblemName.chk#
# to ProblemName.plt.# and ProblemName.chk.# for
# compatability with ParaView's AMReX reader.

# Based on Workflow/AMReX/Utilities/Files.py

# Example (directory structure): ScratchDir/SingleRunDir/ 
# Exmaple (command): python3 ScratchDir SingleRunDir ProblemName

ScratchDir = argv[1] # The path to the directory containing one or more runs.
RunDir = argv[2] # The directory (relative to ScratchDir) containing a single run.
PlotBaseName = argv[3] + '.plt'
CheckBaseName = argv[3] + '.chk'
fileArray = []

NP = len( PlotBaseName )
NC = len( CheckBaseName )

chdir( ScratchDir + '/' + RunDir )

fileArray \
  = listdir( './'  )

for iFile in range( len( fileArray ) ):

    oFile = fileArray[iFile]

    if ( oFile[0:NP] == PlotBaseName and oFile[NP].isdigit() ):
       
        NewName = argv[3] + '.plt.' + str( oFile[NP:] )

        rename( './' + oFile, './' + NewName )
    
    if ( oFile[0:NC] == CheckBaseName and oFile[NC].isdigit() ):

        NewName = argv[3] + '.chk.' + str( oFile[NC:] )

        rename( './' + oFile, './' + NewName )

