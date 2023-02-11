#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from UtilitiesModule import GetData, GetNorm

"""

Default usage, plots last Plotfile in PlotfileDirectory:

  $ python3 plotFieldsAMReX.py

Alernate usage, plot specific file in PlotfileDirectory:

  $ python3 plotFieldsAMReX.py 10

  will plot the *00000010 Plotfile

"""

#### ========== User Input ==========

# Specify name of problem
ProblemName = 'Advection1D'

# Specify title of figure
FigTitle = ProblemName

# Specify directory containing amrex Plotfiles
PlotfileDirectory = 'thornado/SandBox/AMReX/'
#PlotfileDirectory \
#  = '/home/kkadoogan/Work/Codes/\
#thornado/SandBox/AMReX/Euler_Relativistic_IDEAL_MR/'

# Specify plot file base name
PlotfileBaseName = ProblemName + '.plt'

# Specify field to plot
Field = 'PF_D'

# Specify to plot in log-scale
UseLogScale_X  = False
UseLogScale_Y  = False
UseLogScale_2D = False

# Specify whether or not to use physical units
UsePhysicalUnits = False

# Specify coordinate system (currently supports 'cartesian' and 'spherical')
CoordinateSystem = 'cartesian'

# Max level of refinement to plot (-1 plots leaf elements)
MaxLevel = -1

# Write extra info to screen
Verbose = True

# Use custom limts for y-axis (1D) or colorbar (2D)
UseCustomLimits = False
vmin = 0.0
vmax = 2.0

# Save figure (True) or plot figure (False)
SaveFig = False

# Specify colormap (2D only)
cmap = 'viridis'

#### ====== End of User Input =======

ID      = '{:s}_{:s}'.format( ProblemName, Field )
FigName = 'fig.{:s}.png'.format( ID )

# Append "/" to PlotfileDirectory, if not present
if( not PlotfileDirectory[-1] == '/' ): PlotfileDirectory += '/'

TimeUnit   = '[]'
LengthUnit = '[]'
if( UsePhysicalUnits ):

    TimeUnit   = '[ms]'
    LengthUnit = '[km]'

Data, DataUnit, X1_C, X2_C, X3_C, dX1, dX2, dX3, xL, xH, nX, Time \
  = GetData( PlotfileDirectory, PlotfileBaseName, Field, \
             CoordinateSystem, UsePhysicalUnits, argv = argv, \
             MaxLevel = MaxLevel, \
             ReturnTime = True, ReturnMesh = True, Verbose = True )

nDims = 1
if nX[1] > 1: nDims += 1
if nX[2] > 1: nDims += 1

# Re-define nX
nX[0] = X1_C.shape[0]
nX[1] = X1_C.shape[1]
nX[2] = X1_C.shape[2]

if nDims == 1:

    X1_C = np.copy( X1_C[:,0,0] )
    Data = np.copy( Data[:,0,0] )

    fig, ax = plt.subplots( 1, 1 )

    ax.plot( X1_C, Data, 'k.' )
    if( UseLogScale_X ): ax.set_xscale( 'log' )
    if( UseLogScale_Y ): ax.set_yscale( 'log' )
    if( UseCustomLimits ): ax.set_ylim( vmin, vmax )
    ax.set_xlim( xL[0], xH[0] )
    ax.set_xlabel( 'X1_C' + ' ' + LengthUnit )
    ax.set_ylabel( Field + ' ' + DataUnit )
    ax.grid()

elif( nDims == 2 ):

    X1_C = np.copy( X1_C[:,:,0] )
    X2_C = np.copy( X2_C[:,:,0] )
    X3_C = np.copy( X3_C[:,:,0] )
    dX1  = np.copy( dX1 [:,:,0] )
    dX2  = np.copy( dX2 [:,:,0] )
    dX3  = np.copy( dX3 [:,:,0] )
    Data = np.copy( Data[:,:,0] )

    # pcolormesh wants the corners of the elements
    X1c = np.empty( (nX[0]+1,nX[1]+1), np.float64 )
    X2c = np.empty( (nX[0]+1,nX[1]+1), np.float64 )
    for iX1 in range( nX[0] ):
        for iX2 in range( nX[1] ):
            X1c[iX1,iX2] = X1_C[iX1,iX2] - 0.5 * dX1[iX1,iX2]
            X2c[iX1,iX2] = X2_C[iX1,iX2] - 0.5 * dX2[iX1,iX2]
            if   iX2 == nX[1]-1 and iX1 == nX[0]-1:
                X1c[iX1,iX2+1  ] = X1_C[iX1,iX2] - 0.5 * dX1[iX1,iX2]
                X2c[iX1,iX2+1  ] = X2_C[iX1,iX2] + 0.5 * dX2[iX1,iX2]
                X1c[iX1+1,iX2  ] = X1_C[iX1,iX2] + 0.5 * dX1[iX1,iX2]
                X2c[iX1+1,iX2  ] = X2_C[iX1,iX2] - 0.5 * dX2[iX1,iX2]
                X1c[iX1+1,iX2+1] = X1_C[iX1,iX2] + 0.5 * dX1[iX1,iX2]
                X2c[iX1+1,iX2+1] = X2_C[iX1,iX2] + 0.5 * dX2[iX1,iX2]
            elif iX2 == nX[1]-1:
                X1c[iX1,iX2+1] = X1_C[iX1,iX2] - 0.5 * dX1[iX1,iX2]
                X2c[iX1,iX2+1] = X2_C[iX1,iX2] + 0.5 * dX2[iX1,iX2]
            elif iX1 == nX[0]-1:
                X1c[iX1+1,iX2] = X1_C[iX1,iX2] + 0.5 * dX1[iX1,iX2]
                X2c[iX1+1,iX2] = X2_C[iX1,iX2] - 0.5 * dX2[iX1,iX2]

    '''
    # To make lineout plot
    # From: https://yt-project.org/doc/visualizing/
    #       manual_plotting.html#line-plots

    oray = ds.ortho_ray( axis = 0, coords = (0,0) )

    plt.plot( X1_C[:,0], oray[Field] )

    if( UseLogScale_X ): plt.xscale( 'log' )
    if( UseLogScale_Y ): plt.yscale( 'log' )
    plt.show()
    exit()
    '''

    if not UseCustomLimits:

        vmin = Data.min()
        vmax = Data.max()

    Norm = GetNorm( UseLogScale_2D, Data, vmin = vmin, vmax = vmax )

    if CoordinateSystem == 'cartesian':

        fig, ax = plt.subplots( 1, 1 )

        im = ax.pcolormesh( X1c, X2c, Data, \
                            cmap = cmap, \
                            norm = Norm, \
                            shading = 'flat' )

        ax.set_xlim( xL[0], xH[0] )
        ax.set_ylim( xL[1], xH[1] )

        ax.set_xlabel( 'X1_C' + ' ' + LengthUnit )
        ax.set_ylabel( 'X2_C' + ' ' + LengthUnit )

        ax.text( 0.4, 0.9, 'Time = {:.2e} {:}'.format \
                 ( Time, TimeUnit ), \
                   transform = ax.transAxes )

        cbar = fig.colorbar( im )
        cbar.set_label( Field + DataUnit )

    elif CoordinateSystem == 'spherical':

        fig = plt.figure()

        ax = fig.add_subplot( 111, polar = True )

        im = ax.pcolormesh( X2c, X1c, Data, \
                            cmap = cmap, \
                            norm = Norm, \
                            shading = 'flat' )

        ax.set_thetamin( 0.0  )
        ax.set_thetamax( 180.0)
        ax.set_theta_direction( -1 )
        ax.set_theta_zero_location( 'W' ) # z-axis horizontal

ax.set_title( FigTitle )

if SaveFig:

    plt.savefig( FigName, dpi = 300 )

else:

    plt.show()

plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
