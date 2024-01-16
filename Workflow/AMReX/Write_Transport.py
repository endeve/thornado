#!/usr/local/bin/python3

import yt
import numpy as np
import subprocess
from os import listdir
from os.path import isfile
from sys import argv, exit
import matplotlib.pyplot as plt
import math



nsp=2048
nEn=32



"""
Default use (python3 PlotAMReX_yt.py), plots last plot-file in ProblemDirectory

Can also specify certain plot-file:  python3 PlotAMReX_yt.py thornado_00000010
"""

yt.funcs.mylog.setLevel(0) # Suppress initial yt output to screen

# --- Get user's HOME directory ---
HOME = subprocess.check_output( ["echo $HOME"], shell = True)
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# --- Get user's THORNADO_DIR directory ---
THORNADO_DIR = subprocess.check_output( ["echo $THORNADO_DIR"], shell = True)
THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

#### ========== User Input ==========

# Specify directory containing plotfiles
ProblemDirectory = THORNADO_DIR + 'SandBox/AMReX/TwoMoment_Test/'

# Specify name of problem (only used for name of output file(s))
ProblemName = 'KHI'

# Specify plot file base name
PlotFileBaseName = 'thornado'

# Specify field to plot
VariableToPlot = 'PR_D'

# Specify to plot in log-scale
UseLogScale = False

# Specify whether or not to use physical units
UsePhysicalUnits = True

# Specify coordinate system (currently supports 'cartesian' and 'spherical' )
CoordinateSystem = 'cartesian'

# Specify aspect ratio (relativistic KHI needs aspect = 0.5)
aspect = 1.0

# Specify colormap
cmap = 'jet'

# Specify to write out data to a file
WriteOut = 1
# Specify to read in data from a file
ReadIn = 1

#init or final
time = 1
#number of cells
cells = "8"

# Specify whether or not to make a movie
MakeMovie, DataFileName, TimeFileName = False, 'MovieData.dat', 'MovieTime.dat'

#### ================================

if( len( argv ) > 1 ):
    File = argv[1]
else:
    # Get last plotfile in directory
    FileArray \
      = np.sort(np.array( [ file for file in listdir( ProblemDirectory ) ] ))
    FileList = []
    for iFile in range( FileArray.shape[0] ):
        sFile = FileArray[iFile]
        if( sFile[0:len(PlotFileBaseName)+1] == PlotFileBaseName + '_' \
              and sFile[len(PlotFileBaseName)+1].isdigit() ):
            FileList.append( sFile )
    FileArray = np.array( FileList )
    File = FileArray[-1]

# Remove "/" at end of filename, if present
if ( File[-1] == '/' ): File = File[:-1]

ds = yt.load( '{:}'.format( ProblemDirectory + File ) )

print( 'Reading from file: {:}'.format( File ) )
MaxLevel = ds.index.max_level
Time     = ds.current_time
nX       = ds.domain_dimensions
xL       = ds.domain_left_edge
xH       = ds.domain_right_edge

# Get dimensionality of problem
if  ( nX[1] == 1 and nX[2] == 1 ):
    nDims = 1
elif( nX[2] == 1 ):
    nDims = 2
else:
    nDims = 3

CoveringGrid \
  = ds.covering_grid \
      ( level           = MaxLevel, \
        left_edge       = xL, \
        dims            = nX * 2**MaxLevel, \
        num_ghost_zones = nX[0] )

# XXX.to_ndarray() strips array of yt units

D=np.zeros([nEn,nsp],dtype=np.float64)
for i in range(0,nEn):
    n=i+1
    if (n<10):
        name="00"+str(n)
    elif (n>=10 and n<100):
        name="0"+str(n)
    else:
        name=str(n)
    Data = CoveringGrid['PR_D_'+name+'_001' ].to_ndarray()
    for j in range(0,nsp):
        D[i][j]=Data[j];
        
np.savetxt('D.txt', D, delimiter=',');


I1=np.zeros([nEn,nsp],dtype=np.float64)
for i in range(0,nEn):
    n=i+1
    if (n<10):
        name="00"+str(n)
    elif (n>=10 and n<100):
        name="0"+str(n)
    else:
        name=str(n)
    Data = CoveringGrid['PR_I1_'+name+'_001' ].to_ndarray()
    for j in range(0,nsp):
        I1[i][j]=Data[j];

np.savetxt('I1.txt', I1, delimiter=',');





N=np.zeros([nEn,nsp],dtype=np.float64)
for i in range(0,nEn):
    n=i+1
    if (n<10):
        name="00"+str(n)
    elif (n>=10 and n<100):
        name="0"+str(n)
    else:
        name=str(n)
    Data = CoveringGrid['CR_N_'+name+'_001' ].to_ndarray()
    for j in range(0,nsp):
        N[i][j]=Data[j];
        
np.savetxt('N.txt', N, delimiter=',');


G1=np.zeros([nEn,nsp],dtype=np.float64)
for i in range(0,nEn):
    n=i+1
    if (n<10):
        name="00"+str(n)
    elif (n>=10 and n<100):
        name="0"+str(n)
    else:
        name=str(n)
    Data = CoveringGrid['CR_G1_'+name+'_001' ].to_ndarray()
    for j in range(0,nsp):
        G1[i][j]=Data[j];

np.savetxt('G1.txt', G1, delimiter=',');






X1=np.zeros(nsp,dtype=np.float64)
V1=np.zeros(nsp,dtype=np.float64)

Data = CoveringGrid['X1_C' ].to_ndarray()
for k in range(0,nsp):
    X1[k]=Data[k]
    
Data = CoveringGrid['PF_V1' ].to_ndarray()
for k in range(0,nsp):
    V1[k]=Data[k]
np.savetxt('V1.txt', V1, delimiter=',');


np.savetxt('X1.txt', X1, delimiter=',');





RMS2=np.zeros(nsp,dtype=np.float64)
Data = CoveringGrid['GR_RMS_001' ].to_ndarray()
for k in range(0,nsp):
    RMS2[k]=Data[k]
np.savetxt('RMS.txt', RMS2, delimiter=',');


FF=np.zeros(nsp,dtype=np.float64)
Data = CoveringGrid['GR_F_001' ].to_ndarray()
for k in range(0,nsp):
    FF[k]=Data[k]
np.savetxt('FF.txt', FF, delimiter=',');


Chi=np.zeros(nsp,dtype=np.float64)
Data = CoveringGrid['GR_K_001' ].to_ndarray()
for k in range(0,nsp):
    Chi[k]=Data[k]
np.savetxt('Chi.txt', FF, delimiter=',');

IntI1=np.zeros(nsp,dtype=np.float64)
Data = CoveringGrid['GR_I1_001' ].to_ndarray()
for k in range(0,nsp):
    IntI1[k]=Data[k]
np.savetxt('IntI1.txt', IntI1, delimiter=',');

IntegratedD=np.zeros(nsp,dtype=np.float64)
Data = CoveringGrid['GR_D_001' ].to_ndarray()
for k in range(0,nsp):
    IntegratedD[k]=Data[k]
np.savetxt('IntD.txt', IntegratedD, delimiter=',');


IntH1=np.zeros(nsp,dtype=np.float64)
Data = CoveringGrid['GR_H1_001' ].to_ndarray()
for k in range(0,nsp):
    IntH1[k]=Data[k]
np.savetxt('IntH1.txt', IntH1, delimiter=',');

IntJ=np.zeros(nsp,dtype=np.float64)
Data = CoveringGrid['GR_J_001' ].to_ndarray()
for k in range(0,nsp):
    IntJ[k]=Data[k]
np.savetxt('IntJ.txt', IntH1, delimiter=',');



RMSalpha=np.zeros(nsp,dtype=np.float64)

N1=np.zeros(nsp,dtype=np.float64)
F1=np.zeros(nsp,dtype=np.float64)

alpha=np.zeros(nsp,dtype=np.float64)
Data = CoveringGrid['GF_Alpha' ].to_ndarray()
for k in range(0,nsp):
    alpha[k]=Data[k]
np.savetxt('alpha.txt', alpha, delimiter=',');



Psi=np.zeros(nsp,dtype=np.float64)
Data = CoveringGrid['GF_Psi' ].to_ndarray()
for k in range(0,nsp):
    Psi[k]=Data[k]
np.savetxt('Psi.txt', Psi, delimiter=',');


SqrtGm=np.zeros(nsp,dtype=np.float64)
Data = CoveringGrid['GF_SqrtGm' ].to_ndarray()
for k in range(0,nsp):
    SqrtGm[k]=Data[k]
np.savetxt('SqrtGm.txt', SqrtGm, delimiter=',');


for i in range(0,nsp):
    RMSalpha[i]=RMS2[i]*alpha[i]
    N1[i]=SqrtGm[i]*alpha[i]*IntI1[i]
    F1[i]=SqrtGm[i]*alpha[i]**2*IntH1[i]
    
np.savetxt('RMSalpha.txt', RMSalpha, delimiter=','); 
np.savetxt('N1.txt', N1, delimiter=','); 
np.savetxt('F1.txt', F1, delimiter=','); 
