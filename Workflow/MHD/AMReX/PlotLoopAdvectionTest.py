#!/usr/bin/env python3

import numpy as np
import subprocess
from sys import argv
import matplotlib.pyplot as plt
from matplotlib import animation

from UtilitiesModule import GetData, GetNorm, ChoosePlotFile
from MakeDataFile import MakeDataFile

plt.rcParams.update({'font.size': 11})

# --- Get user's HOME directory ---

HOME = subprocess.check_output( ["echo $HOME"], shell = True )
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# --- Get user's THORNADO_DIR directory ---

THORNADO_DIR = subprocess.check_output( ["echo $THORNADO_DIR"], shell = True )
THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

#### ========== User Input ==========

# Specify Problem ID for Figure Saving

ProblemID = 'LoopAdvection2D_VxOnly_001Periods_Con2Prim1.00e-08_1stOrder_n512x512x001_Optimize_Parallel'

# Specify problem directory.

ProblemDir = THORNADO_DIR + 'SandBox/AMReX/MHD_Relativistic_IDEAL/LoopAdvection2D_VxOnly_001Periods_Con2Prim1.00e-08_1stOrder_n512x512x001_Optimize_Parallel/Plotfiles'

# Specify title of figure

FigTitle = 'Loop Advection, 1st Order, 512x512x001 Elements \n' \
           + 'No Divergence Cleaning \n' \
           + 't = 2.00e+00'

# Specify plot file base name

PlotFileBaseName = 'plt'

# Specify field to plot

Field = 'CM_B1'

# Specify colormap

Verbose = True

# Specify saving of figure, figure directory, and figure name.

SaveFig = True

FigDir = ProblemDir + '../' + 'Images/'

FigName = './' + ProblemID + '_' + Field

# Specify use of physical units.

UsePhysicalUnits = True

# Specify Physical Dimensionality of Loop Advection Problem

ProblemDim = '2D'

# Specify if a movie is to be made.

MakeMovie = False

# Make 1D Plot

OneDPlot = False

#### ====== End of User Input =======

def UpdateFrame(t):

    im.set_array()

CoordinateSystem = 'cartesian'

if( not ProblemDir[-1] == '/' ): ProblemDir += '/'

if MakeMovie:

   File, FileArray = ChoosePlotFile( ProblemDir, PlotFileBaseName, \
                                     argv, ReturnFileArray = True, \
                                     Verbose = True )


else:

    Data, DataUnit, X1, X2, X3, xL, xU, nX, Time \
      = GetData( ProblemDir, PlotFileBaseName, Field, \
                 CoordinateSystem, UsePhysicalUnits, argv = argv, \
                 ReturnTime = True, ReturnMesh = True, Verbose = Verbose )

    nDims = 1
    if nX[1] > 1: nDims += 1
    if nX[2] > 1: nDims += 1

    vmin =  None
    vmax =  None

    if( ProblemDim == '3D' ):
        plt.imshow(Data, cmap = 'Greys' )
    else:

        if( Field == 'CM_B1' ):

            vmin   = -10**-3
            vmax   =  10**-3
            clabel = 'Eulerian Magnet'

        elif( Field == 'EulerianMagneticPressure' or Field == 'MagneticPressure' ):

            print(np.max(Data))

            vmin   = 10**-7
            vmax   = 10**-6
            clabel = 'Magnetic Field Pressure'

        if( OneDPlot ):

            plt.plot(X2[64,:], Data[64,:], linewidth = 2)

            plt.xlabel('y @ x = 0.00390625')
            plt.ylabel('Eulerian B2', loc = 'center')
            plt.title(FigTitle + ', t = ' + str(Time))

            plt.xlim([-0.51, 0.51])
            plt.ylim([-10**-3, 10**-3])

            plt.tight_layout()

            FigName = FigName + '_Halfway'

        else:

            plt.figure(figsize = (7.2, 5.4))

            plt.pcolormesh(X2, X1, Data.transpose(), cmap = 'Greys', \
                           vmin = vmin, vmax = vmax)

            plt.xlabel('x')
            plt.ylabel('y')
            plt.xlim(-0.5, 0.5)
            plt.ylim(-0.5, 0.5)
            cbar = plt.colorbar(format = '%.1e')
            cbar.set_label('Eulerian B1')
            plt.title(FigTitle)

    if SaveFig:

        plt.savefig(FigName + '.png', format = 'png')

    plt.show()
