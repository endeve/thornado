# Code developed off of PlotFieldsAMReX.py

#!/usr/bin/env python3

import numpy as np
import subprocess
from sys import argv
import matplotlib.pyplot as plt
from matplotlib import animation

from UtilitiesModule import GetData, GetNorm
from MakeDataFile import MakeDataFile

# --- Get user's HOME directory ---

HOME = subprocess.check_output( ["echo $HOME"], shell = True )
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# --- Get user's THORNADO_DIR directory ---

THORNADO_DIR = subprocess.check_output( ["echo $THORNADO_DIR"], shell = True )
THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

# Specify problem dimensionality.

ProblemDim = 3

# Specify Problem ID for Figure Saving

ProblemID = 'HydroAdvection3D_001Periods_Prim2Con1.00e-16_2ndOrder_n128x001x001_Optimize_Parallel'

# Specify problem directory.

HydroDirList = \
[ 'SandBox/AMReX/Euler_Relativistic_IDEAL/HydroAdvection3D_001Periods_Prim2Con1.00e-16_1stOrder_n002x002x128_Optimize_Parallel/Plotfiles', \
  'SandBox/AMReX/Euler_Relativistic_IDEAL/HydroAdvection3D_001Periods_Prim2Con1.00e-16_2ndOrder_n002x002x128_Optimize_Parallel/Plotfiles', \
  'SandBox/AMReX/Euler_Relativistic_IDEAL/HydroAdvection3D_001Periods_Prim2Con1.00e-16_3rdOrder_n002x002x128_Optimize_Parallel/Plotfiles', \
  'SandBox/AMReX/Euler_Relativistic_IDEAL/HydroAdvection3D_001Periods_Prim2Con1.00e-08_3rdOrder_n002x002x128_Optimize_Parallel/Plotfiles' ]

MHDDirList = \
[ 'SandBox/AMReX/MHD_Relativistic_IDEAL/HydroAdvection3D_001Periods_Prim2Con1.00e-16_1stOrder_n002x002x128_Optimize_Parallel/Plotfiles',  \
  'SandBox/AMReX/MHD_Relativistic_IDEAL/HydroAdvection3D_001Periods_Prim2Con1.00e-16_2ndOrder_n002x002x128_Optimize_Parallel/Plotfiles', \
  'SandBox/AMReX/MHD_Relativistic_IDEAL/HydroAdvection3D_001Periods_Prim2Con1.00e-16_3rdOrder_n002x002x128_Optimize_Parallel/Plotfiles', \
  'SandBox/AMReX/MHD_Relativistic_IDEAL/HydroAdvection3D_001Periods_Prim2Con1.00e-08_3rdOrder_n002x002x128_Optimize_Parallel/Plotfiles' ]

# Specify title of figure
FigTitle = '3D Hydro Sine Wave Advection Comparison (Hydro Code vs. MHD Code) \n' \
           + '1 Period, 002x002x128 Elements \n'

# Specify plotfile base name

PlotFileBaseName = 'plt'

# Specify fields to compare.

HydroField = 'CF_D'

MHDField   = 'CM_D'

# Specify verbosity

Verbose = True

# Specify saving of figure.

SaveFig = True

# Specify use of physical units.

UsePhysicalUnits = True

# Specify if a movie is to be made.

MakeMovie = False

# Specify linestyles and colors.

Linestyles = ['-', '-', '-', '--']

Colors = ['b', 'g', 'r', 'r']

#### ====== End of User Input =======

CoordinateSystem = 'cartesian'

ID           = '{:s}_{:s}'.format( ProblemID, 'P_D' )
FigName      = 'fig.{:s}.png'.format( ID )

for i, HydroDir in enumerate(HydroDirList):
    HydroDirList[i] = THORNADO_DIR + HydroDir
    if( not HydroDir[-1] == '/' ): HydroDirList[i] += '/'

for i, MHDDir in enumerate(MHDDirList):
    MHDDirList[i] = THORNADO_DIR + MHDDir
    if( not MHDDir[-1] == '/' ): MHDDirList[i] += '/'

HydroDataList      = []
HydroDataUnitList  = []

MHDDataList     = []
MHDDataUnitList = []

ComparisonDataList = []

JunkData, JunkDataUnit, X1, X2, X3, xL, xU, nX, Time \
= GetData(HydroDirList[0], PlotFileBaseName, HydroField, CoordinateSystem, \
          UsePhysicalUnits, argv = argv, ReturnTime = True, ReturnMesh = True, Verbose = Verbose)

for i, HydroDir, MHDDir in zip(range(0, len(HydroDirList)), HydroDirList, MHDDirList):

    HydroDataList.append(GetData(HydroDir, PlotFileBaseName, HydroField, \
                         CoordinateSystem, UsePhysicalUnits, argv = argv, \
                         ReturnTime = True, ReturnMesh = False, Verbose = Verbose)[0])

    MHDDataList.append(GetData(MHDDir, PlotFileBaseName, MHDField, \
                       CoordinateSystem, UsePhysicalUnits, argv = argv, \
                       ReturnTime = True, ReturnMesh = False, Verbose = Verbose)[0])

    ComparisonDataList.append( np.abs( HydroDataList[i] - MHDDataList[i] ) / HydroDataList[i] )

plt.figure(2)

for i, ComparisonData in enumerate(ComparisonDataList):

   if( ProblemDim == 1 ):

       plt.plot(X1, ComparisonData)

   elif( ProblemDim == 2 ):

       plt.plot(X2, np.max(ComparisonData, axis = 0))

   elif( ProblemDim == 3 ):

       plt.semilogy(X3[0,:,0], np.max( np.max( ComparisonData, axis = 0), axis = 0), linewidth = 2.5, linestyle = Linestyles[i], color = Colors[i])
       plt.xlim([0, 1.0])

plt.xlabel('X3')
plt.ylabel('Absolute Relative Error (max in X1 and X2 directions)')
plt.legend(['Order 1 - P2C Tols. = $10^{-16}$', 'Order 2 - P2C Tols. = $10^{-16}$', \
            'Order 3 - P2C Tols. = $10^{-16}$', 'Order 3 - P2C Tols. = $10^{-8}$'], \
            loc = 'center left', bbox_to_anchor=(0.6, 0.65), fontsize = 9)

plt.savefig('HydroAdvection3DComparison_001Periods_n002x002x128_Optimize_Parallel.png')

plt.close()
