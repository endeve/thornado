#!/usr/bin/env python3

import numpy as np
import subprocess
from sys import argv
import matplotlib.pyplot as plt

from UtilitiesModule import GetData, GetNorm
from MakeDataFile import MakeDataFile

# --- Get user's HOME directory ---

HOME = subprocess.check_output( ["echo $HOME"], shell = True )
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# --- Get user's THORNADO_DIR directory ---

THORNADO_DIR = subprocess.check_output( ["echo $THORNADO_DIR"], shell = True )
THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

# Specify problem directories.

DirList = \
[
  'SandBox/AMReX/MHD_Relativistic_IDEAL/NewCPAlfvenX3ConvergenceTest/0008', \
  'SandBox/AMReX/MHD_Relativistic_IDEAL/NewCPAlfvenX3ConvergenceTest/0016', \
  'SandBox/AMReX/MHD_Relativistic_IDEAL/NewCPAlfvenX3ConvergenceTest/0032', \
  'SandBox/AMReX/MHD_Relativistic_IDEAL/NewCPAlfvenX3ConvergenceTest/0064', \
  'SandBox/AMReX/MHD_Relativistic_IDEAL/NewCPAlfvenX3ConvergenceTest/0128', \
  'SandBox/AMReX/MHD_Relativistic_IDEAL/NewCPAlfvenX3ConvergenceTest/0256', \
  'SandBox/AMReX/MHD_Relativistic_IDEAL/NewCPAlfvenX3ConvergenceTest/0512', \
  'SandBox/AMReX/MHD_Relativistic_IDEAL/NewCPAlfvenX3ConvergenceTest/1024', \
  'SandBox/AMReX/MHD_Relativistic_IDEAL/NewCPAlfvenX3ConvergenceTest/2048'
]

#DirList = \
#[
#  'SandBox/AMReX/MHD_Relativistic_IDEAL/CPAlfven3D_001Periods_Con2Prim1.00e-08_1stOrder_n002x002x008_Optimize_Parallel/Plotfiles', \
#  'SandBox/AMReX/MHD_Relativistic_IDEAL/CPAlfven3D_001Periods_Con2Prim1.00e-08_1stOrder_n002x002x016_Optimize_Parallel/Plotfiles', \
#  'SandBox/AMReX/MHD_Relativistic_IDEAL/CPAlfven3D_001Periods_Con2Prim1.00e-08_1stOrder_n002x002x032_Optimize_Parallel/Plotfiles', \
#  'SandBox/AMReX/MHD_Relativistic_IDEAL/CPAlfven3D_001Periods_Con2Prim1.00e-08_1stOrder_n002x002x064_Optimize_Parallel/Plotfiles', \
#  'SandBox/AMReX/MHD_Relativistic_IDEAL/CPAlfven3D_001Periods_Con2Prim1.00e-08_1stOrder_n002x002x128_Optimize_Parallel/Plotfiles'
#]

# Specify specific plot file numbers.

FileNumberList = [ '', '', '', '', '', '', '',  '', '' ]

# Specify plotfile base name

PlotFileBaseName = 'plt'

# Specify field to calculate error

Field = 'CM_B2'

# Specify verbosity

Verbose = True

# Specify use of physical units

UsePhysicalUnits = False

# Specify Order

Order = 1

#### ====== End of User Input =======

CoordinateSystem = 'cartesian'

#ID           = '{:s}_{:s}'.format( ProblemID, 'P_D' )
#FigName      = 'fig.{:s}.png'.format( ID )

for i, Dir in enumerate(DirList):
    DirList[i] = THORNADO_DIR + Dir
    if( not Dir[-1] == '/' ): DirList[i] += '/'

nDOFXList    = []

L1ErrorList   = []
L2ErrorList   = []
LInfErrorList = []

E = []

for i, Dir in zip(range(0, len(DirList)), DirList):

    InitialData, DataUnit, X1, X2, X3, xL, xU, nX, Time \
    = GetData(Dir, PlotFileBaseName, Field, CoordinateSystem, \
              UsePhysicalUnits, argv = [argv, '0'], ReturnTime = True, ReturnMesh = True, Verbose = Verbose)

    Normalization = np.max( np.abs( InitialData ) )

    print(Normalization)

    if FileNumberList[i] != '':

        FinalData, DataUnit, X1, X2, X3, xL, xU, nX, Time \
        = GetData(Dir, PlotFileBaseName, Field, CoordinateSystem, \
                  UsePhysicalUnits, argv = [argv, FileNumberList[i]], ReturnTime = True, ReturnMesh = True, Verbose = Verbose)

    else:

        FinalData, DataUnit, X1, X2, X3, xL, xU, nX, Time \
        = GetData(Dir, PlotFileBaseName, Field, CoordinateSystem, \
                  UsePhysicalUnits, argv = argv, ReturnTime = True, ReturnMesh = True, Verbose = Verbose)

    L1ErrorList.append( np.sum( np.abs( FinalData - InitialData ) / Normalization ) / ( nX[0] * nX[1] * nX[2] ) )
    L2ErrorList.append( np.sqrt( np.sum( ( FinalData - InitialData )**2 ) / ( nX[0] * nX[1] * nX[2] ) ) / Normalization )
    LInfErrorList.append( np.max( np.abs( FinalData - InitialData ) ) / Normalization )

    nDOFXList.append(nX[0] * nX[1] * nX[2])

for i in range(0, len(nDOFXList)):

    E.append( ( LInfErrorList[-1] * nDOFXList[-1]**Order ) / nDOFXList[i]**Order )

plt.loglog(nDOFXList, LInfErrorList, marker = '+', color = 'g', linewidth = 2)
plt.loglog(nDOFXList, L2ErrorList,   marker = 'x', color = 'b', linewidth = 2)
plt.loglog(nDOFXList, L1ErrorList,   marker = 's', color = 'r', linewidth = 2)

plt.loglog(nDOFXList, E, color = 'k', linewidth = 2)

plt.rcParams.update({'font.size': 11, 'lines.linewidth': 2})
plt.xlabel('# of Elements (2 x 2 x X)')
plt.ylabel('Normalized Error (Eulerian B2)')
plt.legend([r'$L_{Inf}$', r'$L_{2}$', r'$L_{1}$', r'Order ' + str(Order) + ' Rate'])
plt.title('Circularly Polarized Alfvén Wave (Z-Axis) \n' + '1st Order, 002x002xX Elements, 1 Period')
plt.savefig('Convergence_CPAlfvenX3_0001Periods_1stOrder_Debug_Parallel_CB2.png')

plt.show()
