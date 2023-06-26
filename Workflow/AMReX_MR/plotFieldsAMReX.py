#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from UtilitiesModule import GetData, GetNorm, MapCenterToCorners

"""

Default usage, plots last Plotfile in PlotfileDirectory:

  $ python3 plotFieldsAMReX.py

Alernate usage, plot specific file in PlotfileDirectory:

  $ python3 plotFieldsAMReX.py 10

  will plot the *00000010 Plotfile

"""

#### ========== User Input ==========

# Specify name of problem
ProblemName = 'Advection2D'

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

polar = False
if CoordinateSystem == 'spherical':
    polar = True

ID      = '{:s}_{:s}'.format( ProblemName, Field )
FigName = 'fig.{:s}.png'.format( ID )

# Append "/" to PlotfileDirectory, if not present
if not PlotfileDirectory[-1] == '/': PlotfileDirectory += '/'

TimeUnits = ''
X1Units   = ''
X2Units   = ''
if UsePhysicalUnits:
    TimeUnits = 'ms'
    if   CoordinateSystem == 'cartesian':
        X1Units = 'km'
        X2Units = 'km'
    elif CoordinateSystem == 'spherical':
        X1Units = 'km'
        X2Units = 'rad'

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

    dX1  = np.copy( dX1 [:,0,0] )
    X1_C = np.copy( X1_C[:,0,0] )
    Data = np.copy( Data[:,0,0] )

    fig, ax = plt.subplots( 1, 1 )

    ax.plot( X1_C, Data, 'k.' )
    if UseLogScale_X:
        ax.set_xscale( 'log' )
        xL = [ max( xL[0], 0.0 + 0.25 * dX1[0] ), 0 ]
    if UseLogScale_Y: ax.set_yscale( 'log' )
    if UseCustomLimits: ax.set_ylim( vmin, vmax )
    ax.set_xlim( xL[0], xH[0] )
    ax.set_xlabel \
      ( r'$x^{{1}}\ \left[\mathrm{{{:}}}\right]$'.format( X1Units ), \
        fontsize = 15 )
    ax.set_ylabel( Field + ' ' + DataUnit )
    ax.grid()

elif nDims == 2:

    X1_C = np.copy( X1_C[:,:,0] )
    X2_C = np.copy( X2_C[:,:,0] )
    X3_C = np.copy( X3_C[:,:,0] )
    dX1  = np.copy( dX1 [:,:,0] )
    dX2  = np.copy( dX2 [:,:,0] )
    dX3  = np.copy( dX3 [:,:,0] )
    Data = np.copy( Data[:,:,0] )

    '''
    # Line-out plots
    fig, ax = plt.subplots( 1, 1 )

    nX2 = X1_C.shape[1]

    ax.plot( X1_C[:,0      ], Data[:,0     ], \
             label = r'$\theta/\pi={:.3f}$'.format( X2_C[0,0     ] / np.pi ) )
    ax.plot( X1_C[:,nX2//2 ], Data[:,nX2//2], \
              label = r'$\theta/\pi={:.3f}$'.format( X2_C[0,nX2//2] / np.pi ) )
    ax.plot( X1_C[:,-1     ], Data[:,-1    ], \
              label = r'$\theta/\pi={:.3f}$'.format( X2_C[0,-1    ] / np.pi ) )

    if UseLogScale_X: ax.set_xscale( 'log' )
    if UseLogScale_Y: ax.set_yscale( 'log' )
    ax.legend()

    plt.show()
    exit()
    '''

    if not UseCustomLimits:

        vmin = Data.min()
        vmax = Data.max()

    Norm = GetNorm( UseLogScale_2D, Data, vmin = vmin, vmax = vmax )

    fig = plt.figure()
    ax = fig.add_subplot( 111, polar = polar )

    # pcolormesh wants the corners of the elements
    X1c, X2c = MapCenterToCorners( X1_C, X2_C, dX1, dX2 )

    if CoordinateSystem == 'cartesian':

        ax.set_xlim( xL[0], xH[0] )
        ax.set_ylim( xL[1], xH[1] )

        ax.set_xlabel \
          ( r'$x^{{1}}\ \left[\mathrm{{{:}}}\right]$'.format( X1Units ), \
            fontsize = 15 )
        ax.set_ylabel \
          ( r'$x^{{2}}\ \left[\mathrm{{{:}}}\right]$'.format( X2Units ), \
            fontsize = 15 )

        ax.text( 0.4, 0.9, r'$t={:.2e}\ \left[\mathrm{{{:}}}\right]$'.format \
                 ( Time, TimeUnits ), transform = ax.transAxes )

    elif CoordinateSystem == 'spherical':

        ax.grid( False )
        ax.set_thetamin( 0.0  )
        ax.set_thetamax( 180.0)
        ax.set_theta_direction( -1 )
        ax.set_theta_zero_location( 'W' ) # z-axis horizontal

        ax.text( 0.6, 0.9, r'$t={:.2e}\ \left[\mathrm{{{:}}}\right]$'.format \
                 ( Time, TimeUnits ), transform = ax.transAxes )

        # Transpose data for spherical-polar coordinates

        X1c = np.copy( X1c[:,0] )
        X2c = np.copy( X2c[0,:] )

        X1c, X2c = np.meshgrid( X2c, X1c )

    im = ax.pcolormesh( X1c, X2c, Data, \
                        cmap = cmap, \
                        norm = Norm, \
                        shading = 'flat' )

    cbar = fig.colorbar( im )
    cbar.set_label( Field + ' ' + r'$\left[\mathrm{{{:}}}\right]$' \
                                  .format( DataUnit ) )

ax.set_title( r'$\texttt{{{:}}}$'.format( FigTitle ) )

if SaveFig:

    plt.savefig( FigName, dpi = 300 )

else:

    plt.show()

plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
