#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
plt.style.use( 'publication.sty' )

from UtilitiesModule import GetNorm, MapCenterToCorners
from MakeDataFile import MakeDataFile, ReadHeader

"""

Creates a directory with structure as laid out
in MakeDataFile.py and makes a movie from it

Usage:
  $ python3 makeMovie2D.py

"""

#### ========== User Input ==========

# ID to be used for naming purposes
ID = 'Advection2D'

# Directory containing AMReX plotfiles
plotfileDirectory = 'thornado/SandBox/AMReX/'
#plotfileDirectory \
#  = '/home/kkadoogan/Work/Codes/thornado/\
#SandBox/AMReX/Euler_Relativistic_IDEAL_MR/'

# plotfile base name (e.g., Advection2D.plt######## -> Advection2D.plt )
plotfileBaseName = ID + '.plt'

# Field to plot
Field = 'PF_D'

# Plot data in log10-scale?
UseLogScale = False

# Unit system of the data
UsePhysicalUnits = False

# Coordinate system (currently supports 'cartesian' and 'spherical' )
CoordinateSystem = 'cartesian'

# Only use every <plotEvery> plotfile
plotEvery = 1

# Colormap
cmap = 'viridis'

# First and last snapshots and number of snapshots to include in movie
SSi = -1 # -1 -> SSi = 0
SSf = -1 # -1 -> plotfileArray.shape[0] - 1
nSS = -1 # -1 -> plotfileArray.shape[0]

# Max level of refinement to include
MaxLevel = -1 # -1 -> use all levels

Verbose = True

UseCustomLimits = False
vmin = 0.0
vmax = 2.0

MovieRunTime = 10.0 # seconds

#### ====== End of User Input =======

if CoordinateSystem == 'spherical':
    polar = True
else:
    polar = False

DataDirectory = '.{:s}'.format( ID )
MovieName     = 'mov.{:s}_{:s}.mp4'.format( ID, Field )

# Append "/" if not present
if not plotfileDirectory[-1] == '/': plotfileDirectory += '/'
if not DataDirectory    [-1] == '/': DataDirectory     += '/'

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

plotfileArray \
  = MakeDataFile( Field, plotfileDirectory, DataDirectory, \
                  plotfileBaseName, CoordinateSystem, \
                  SSi = SSi, SSf = SSf, nSS = nSS, \
                  UsePhysicalUnits = UsePhysicalUnits, \
                  MaxLevel = MaxLevel, Verbose = Verbose )
plotfileArray = np.copy( plotfileArray[::plotEvery] )

def f(t):

    FileDirectory = DataDirectory + plotfileArray[t] + '/'

    TimeFile = FileDirectory + '{:}.dat'.format( 'Time' )
    X1File   = FileDirectory + '{:}.dat'.format( 'X1' )
    X2File   = FileDirectory + '{:}.dat'.format( 'X2' )
    dX1File  = FileDirectory + '{:}.dat'.format( 'dX1' )
    dX2File  = FileDirectory + '{:}.dat'.format( 'dX2' )
    DataFile = FileDirectory + '{:}.dat'.format( Field )

    DataShape, DataUnits, MinVal, MaxVal = ReadHeader( DataFile )

    Time = np.loadtxt( TimeFile )
    X1_C = np.loadtxt( X1File   ).reshape( np.int64( DataShape ) )
    X2_C = np.loadtxt( X2File   ).reshape( np.int64( DataShape ) )
    dX1  = np.loadtxt( dX1File  ).reshape( np.int64( DataShape ) )
    dX2  = np.loadtxt( dX2File  ).reshape( np.int64( DataShape ) )
    Data = np.loadtxt( DataFile )

    return Data, DataUnits, X1_C, X2_C, dX1, dX2, Time

Data, DataUnits, X1_C, X2_C, dX1, dX2, Time = f(0)

if nSS < 0: nSS = plotfileArray.shape[0]

if not UseCustomLimits:
    vmin = +np.inf
    vmax = -np.inf
    for j in range( nSS ):
        DataFile \
          = DataDirectory + plotfileArray[j] + '/{:}.dat'.format( Field )
        DataShape, DataUnits, MinVal, MaxVal = ReadHeader( DataFile )
        vmin = min( vmin, MinVal )
        vmax = max( vmax, MaxVal )

nX = np.shape(X1_C)

x1L = X1_C[0 ,0 ] - 0.5 * dX1[0 ,0 ]
x1H = X1_C[-1,-1] + 0.5 * dX1[-1,-1]
x2L = X2_C[0 ,0 ] - 0.5 * dX2[0 ,0 ]
x2H = X2_C[-1,-1] + 0.5 * dX2[-1,-1]

fig = plt.figure()
ax  = fig.add_subplot( 111, polar = polar )
ax.set_title( r'$\texttt{{{:}}}$'.format( ID ) )

time_text = ax.text( 0.6, 0.9, '', transform = ax.transAxes )

X1c, X2c = MapCenterToCorners( X1_C, X2_C, dX1, dX2 )

if CoordinateSystem == 'spherical':

    X1c = np.copy( X1c[:,0] )
    X2c = np.copy( X2c[0,:] )

    X1c, X2c = np.meshgrid( X2c, X1c )

    ax.grid( False )

    ax.set_thetamin( 0.0 )
    ax.set_thetamax( 180.0 )
    ax.set_rmin( 0.0 )
    ax.set_rmax( x1H )

    ax.set_theta_zero_location( 'N' ) # z-axis vertical
    ax.set_theta_direction( -1 )

elif CoordinateSystem == 'cartesian':

    ax.set_xlabel \
      ( r'$x^{{1}}\ \left[\mathrm{{{:}}}\right]$'.format( X1Units ), \
        fontsize = 15 )
    ax.set_ylabel \
      ( r'$x^{{2}}\ \left[\mathrm{{{:}}}\right]$'.format( X2Units ), \
        fontsize = 15 )

Norm = GetNorm( UseLogScale, Data, vmin = vmin, vmax = vmax )

im = ax.pcolormesh( X1c, X2c, Data, \
                    cmap = cmap, \
                    norm = Norm, \
                    shading = 'flat' )

time_text = ax.text( 0.4, 0.9, '', transform = ax.transAxes )

cbar = fig.colorbar( im )
cbar.set_label( Field + ' ' + r'$\left[\mathrm{{{:}}}\right]$' \
                              .format( DataUnits[2:-2] ) )

def InitializeFrame():

    Data, DataUnits, X1_C, X2_C, dX1, dX2, Time = f(0)
    im.set_array( Data.flatten() )
    time_text.set_text('')
    ret = ( im, time_text )
    return ret

def UpdateFrame(t):

    print( '    {:}/{:}'.format( t, nSS ) )
    Data, DataUnits, X1_C, X2_C, dX1, dX2, Time = f(t)

    # pcolormesh wants the corners of the elements
    X1c, X2c = MapCenterToCorners( X1_C, X2_C, dX1, dX2 )

    if CoordinateSystem == 'spherical':

        X1c = np.copy( X1c[:,0] )
        X2c = np.copy( X2c[0,:] )

        X1c, X2c = np.meshgrid( X2c, X1c )

    im = ax.pcolormesh( X1c, X2c, Data, \
                        cmap = cmap, \
                        norm = Norm, \
                        shading = 'flat' )

    im.set_array( Data.flatten() )
    time_text.set_text( r'$t={:.3e}\ \left[\mathrm{{{:}}}\right]$' \
                        .format( Time, TimeUnits ) )

    ret = ( im, time_text )

    return ret

# Call the animator
anim \
  = animation.FuncAnimation \
      ( fig, UpdateFrame, init_func = InitializeFrame, \
        frames = nSS, blit = True)

fps = max( 1, nSS // MovieRunTime )

print( '\n  Making Movie' )
print( '  ------------' )

anim.save( MovieName, fps = fps, dpi = 300 )

import os
os.system( 'rm -rf __pycache__ ' )
