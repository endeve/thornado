#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

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
ID = 'KelvinHelmholtz2D'

# Directory containing AMReX Plotfiles
PlotfileDirectory \
  = '/home/kkadoogan/Work/Codes/thornado/SandBox/AMReX/Euler_Relativistic_IDEAL_MR/'

# Plotfile base name (e.g., Advection2D.plt######## -> Advection2D.plt )
PlotfileBaseName = ID + '.plt'

# Field to plot
Field = 'PF_D'

# Plot data in log10-scale?
UseLogScale = False

# Unit system of the data
UsePhysicalUnits = False

# Coordinate system (currently supports 'cartesian' and 'spherical' )
CoordinateSystem = 'cartesian'

# Colormap
cmap = 'jet'

# First and last snapshots and number of snapshots to include in movie
SSi = -1 # -1 -> SSi = 0
SSf = -1 # -1 -> PlotfileArray.shape[0] - 1
nSS = -1 # -1 -> PlotfileArray.shape[0]

# Max level of refinement to include
MaxLevel = -1 # -1 -> use all levels

Verbose = True

UseCustomLimits = False
vmin = 0.0
vmax = 2.0

MovieRunTime = 10.0 # seconds

zAxisVertical = False

#### ====== End of User Input =======

DataDirectory = '.{:s}'.format( ID )
MovieName     = 'mov.{:s}_{:s}.mp4'.format( ID, Field )

# Append "/" if not present
if( not PlotfileDirectory[-1] == '/' ): PlotfileDirectory += '/'
if( not DataDirectory    [-1] == '/' ): DataDirectory     += '/'

TimeUnit = ''
if UsePhysicalUnits: TimeUnit = 'ms'

PlotfileArray \
  = MakeDataFile( Field, PlotfileDirectory, DataDirectory, \
                  PlotfileBaseName, CoordinateSystem, \
                  SSi = SSi, SSf = SSf, nSS = nSS, \
                  UsePhysicalUnits = UsePhysicalUnits, \
                  MaxLevel = MaxLevel, Verbose = Verbose )
def f(t):

    FileDirectory = DataDirectory + PlotfileArray[t] + '/'

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

nX = np.shape(X1_C)

x1L = X1_C[0 ,0 ] - 0.5 * dX1[0 ,0 ]
x1H = X1_C[-1,-1] + 0.5 * dX1[-1,-1]
x2L = X2_C[0 ,0 ] - 0.5 * dX2[0 ,0 ]
x2H = X2_C[-1,-1] + 0.5 * dX2[-1,-1]

if nSS < 0: nSS = PlotfileArray.shape[0]

if not UseCustomLimits:
    vmin = +np.inf
    vmax = -np.inf
    for j in range( nSS ):
        DataFile \
          = DataDirectory + PlotfileArray[j] + '/{:}.dat'.format( Field )
        DataShape, DataUnits, MinVal, MaxVal = ReadHeader( DataFile )
        vmin = min( vmin, MinVal )
        vmax = max( vmax, MaxVal )

Norm = GetNorm( UseLogScale, Data, vmin = vmin, vmax = vmax )

fig = plt.figure()
ax  = fig.add_subplot( 111 )

# pcolormesh wants the corners of the elements
X1c, X2c = MapCenterToCorners( X1_C, X2_C )

im = ax.pcolormesh( X1c, X2c, Data, \
                    cmap = cmap, \
                    norm = Norm, \
                    shading = 'flat' )

ax.set_xlim( x1L, x1H )
ax.set_ylim( x2L, x2H )
#ax.set_aspect( 'auto' )

time_text = plt.text( 0.5, 0.9, '', transform = ax.transAxes )

cbar = fig.colorbar( im )

cbar.set_label( Field + ' ' + DataUnits )

def UpdateFrame(t):

    Data, DataUnits, X1_C, X2_C, dX1, dX2, Time = f(t)

    # pcolormesh wants the corners of the elements
    X1c, X2c = MapCenterToCorners( X1_C, X2_C )

    im = ax.pcolormesh( X1c, X2c, Data, \
                        cmap = cmap, \
                        norm = Norm, \
                        shading = 'flat' )

    im.set_array( Data.flatten() )
    time_text.set_text( 'Time = {:.3e} {:}'.format( Time, TimeUnit ) )

    ret = ( im, time_text )
    return ret

# Call the animator
anim = animation.FuncAnimation \
         ( fig, UpdateFrame, frames = nSS, \
           blit = True)

fps = max( 1, nSS / MovieRunTime )

anim.save( MovieName, fps = fps )

import os
os.system( 'rm -rf __pycache__ ' )
