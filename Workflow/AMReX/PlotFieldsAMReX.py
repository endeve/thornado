#!/usr/bin/env python3

import numpy as np
import subprocess
from sys import argv
import matplotlib.pyplot as plt
from matplotlib import animation

from UtilitiesModule import GetData, GetNorm
from MakeDataFile import MakeDataFile

"""

Default usage, plots last plotfile in DataDirectory:

  $ python3 PlotFieldsAMReX.py

Alernate usage, plot specific file in DataDirectory:

  $ python3 PlotFieldsAMReX.py 10

  will plot the *_00000010 plotfile

"""

# --- Get user's HOME directory ---

HOME = subprocess.check_output( ["echo $HOME"], shell = True )
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# --- Get user's THORNADO_DIR directory ---

THORNADO_DIR = subprocess.check_output( ["echo $THORNADO_DIR"], shell = True )
THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

#### ========== User Input ==========

# Specify name of problem
ProblemName = 'KHI'

# Specify title of figure
FigTitle = ProblemName

# Specify directory containing plotfiles
DataDirectory = THORNADO_DIR + 'SandBox/AMReX/'

# Specify plot file base name
PlotFileBaseName = 'plt'

# Specify field to plot
Field = 'PF_D'

# Specify to plot in log-scale
UseLogScale = False

# Specify whether or not to use physical units
UsePhysicalUnits = False

# Specify coordinate system (currently supports 'cartesian' and 'spherical' )
CoordinateSystem = 'cartesian'

# Specify colormap
cmap = 'jet'

Verbose = True

UseCustomLimits = False
vmin = 0.0
vmax = 2.0

SaveFig = False

#### ====== End of User Input =======

ID           = '{:s}_{:s}'.format( ProblemName, Field )
DataFileName = '.{:s}_MovieData.dat'.format( ID )
FigName      = 'fig.{:s}.png'.format( ID )

# Append "/" to DataDirectory, if not present
if( not DataDirectory[-1] == '/' ): DataDirectory += '/'

TimeUnit = ''
LengthUnit = ''
if( UsePhysicalUnits ):

    TimeUnit = 'ms'
    LengthUnit = 'km'

Data, DataUnit, X1, X2, X3, xL, xU, nX, Time \
  = GetData( DataDirectory, PlotFileBaseName, Field, \
             CoordinateSystem, UsePhysicalUnits, argv = argv, \
             ReturnTime = True, ReturnMesh = True, Verbose = Verbose )

nDims = 1
if nX[1] > 1: nDims += 1
if nX[2] > 1: nDims += 1

assert ( ( nDims >= 1 ) & ( nDims <= 2 ) ), \
       'Invalid nDims: {:d}\nnDims must equal 1 or 2'.format( nDims )

if nDims == 1:

    plt.plot( X1, Data, 'k-' )
    if( UseLogScale ): plt.yscale( 'log' )
    plt.xlim( xL[0], xU[0] )
    plt.xlabel( 'X1' + ' ' + LengthUnit )
    plt.ylabel( Field )
    plt.show()
    plt.close()

elif( nDims == 2 ):

    '''
    # To make lineout plot
    # From: https://yt-project.org/doc/visualizing/
    #       manual_plotting.html#line-plots

    oray = ds.ortho_ray( axis = 0, coords = (0,0) )

    plt.plot( X1, oray[Field] )

    if( UseLogScale ): plt.yscale( 'log' )
    plt.show()
    exit()
    '''

    fig = plt.figure()
    fig.suptitle( FigTitle )

    if not UseCustomLimits:

        vmin = Data.min()
        vmax = Data.max()

    Norm = GetNorm( UseLogScale, Data, vmin = vmin, vmax = vmax )

    if CoordinateSystem == 'spherical':

        ax = fig.add_subplot( 111, polar = True )

        im = ax.pcolormesh( X2, X1, Data, \
                            cmap = cmap, \
                            norm = Norm, \
                            shading = 'nearest' )

        ax.set_thetamin( 0.0  )
        ax.set_thetamax( 180.0)
        ax.set_theta_direction( -1 )
        ax.set_theta_zero_location( 'W' ) # z-axis horizontal

    else:

        ax = fig.add_subplot( 111 )

        im = ax.pcolormesh( X1, X2, Data, \
                            cmap = cmap, \
                            norm = Norm, \
                            shading = 'nearest' )

        ax.set_xlim( xL[0], xU[0] )
        ax.set_ylim( xL[1], xU[1] )

        ax.set_xlabel( 'X1' + ' ' + LengthUnit )
        ax.set_ylabel( 'X2' + ' ' + LengthUnit )

        ax.text( 0.4, 0.9, 'Time = {:.2e} {:}'.format \
                 ( Time, TimeUnit ), \
                   transform = ax.transAxes )

    cbar = fig.colorbar( im )
    cbar.set_label( Field + DataUnit )

    if SaveFig:

        plt.savefig( FigName, dpi = 300 )

    else:

        plt.show()
        plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
