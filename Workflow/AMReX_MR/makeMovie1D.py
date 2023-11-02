#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
plt.style.use( 'publication.sty' )

from MakeDataFile import MakeDataFile, ReadHeader

"""

Creates a directory with structure as laid out
in MakeDataFile.py and makes a movie from it

Usage:
  $ python3 makeMovie1D.py

"""

#### ========== User Input ==========

# ID to be used for naming purposes
ID = 'Advection1D'

# Directory containing AMReX plotfiles
plotfileDirectory = 'thornado/SandBox/AMReX/'
#plotfileDirectory \
#  = '/home/kkadoogan/Work/Codes/thornado/\
#SandBox/AMReX/Euler_Relativistic_IDEAL_MR/'

# plotfile base name (e.g., Advection1D.plt######## -> Advection1D.plt )
plotfileBaseName = ID + '.plt'

# Field to plot
Field = 'PF_D'

# Plot data in log10-scale?
UseLogScale_Y = False
UseLogScale_X = False

# Unit system of the data
UsePhysicalUnits = False

# Coordinate system (currently supports 'cartesian' and 'spherical' )
CoordinateSystem = 'cartesian'

# Only use every <plotEvery> plotfile
plotEvery = 1

# First and last snapshots and number of snapshots to include in movie
SSi = -1 # -1 -> SSi = 0
SSf = -1 # -1 -> plotfileArray.shape[0] - 1
nSS = -1 # -1 -> plotfileArray.shape[0]

# Max level of refinement to include
MaxLevel = -1

# Include initial conditions in movie?
ShowIC = True

PlotMesh = True

Verbose = True

UseCustomLimits = False
vmin = 0.0
vmax = 2.0

MovieRunTime = 10.0 # seconds

#### ====== End of User Input =======

DataDirectory = '.{:s}'.format( ID )
MovieName     = 'mov.{:s}_{:s}.mp4'.format( ID, Field )

# Append "/" if not present
if not plotfileDirectory[-1] == '/': plotfileDirectory += '/'
if not DataDirectory    [-1] == '/': DataDirectory     += '/'

TimeUnits = ''
X1Units   = ''
if UsePhysicalUnits:
    TimeUnits = 'ms'
    X1Units   = 'km'

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
    dX1File  = FileDirectory + '{:}.dat'.format( 'dX1' )
    DataFile = FileDirectory + '{:}.dat'.format( Field )

    DataShape, DataUnits, MinVal, MaxVal = ReadHeader( DataFile )

    Time = np.loadtxt( TimeFile )
    X1_C = np.loadtxt( X1File   )
    dX1  = np.loadtxt( dX1File  )
    Data = np.loadtxt( DataFile )

    return Data, DataUnits, X1_C, dX1, Time

Data0, DataUnits, X1_C0, dX10, Time = f(0)

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

nX = np.shape(X1_C0)

xL = X1_C0[0 ] - 0.5 * dX10[0 ]
xH = X1_C0[-1] + 0.5 * dX10[-1]

fig = plt.figure()
ax  = fig.add_subplot( 111 )
ax.set_title( r'$\texttt{{{:}}}$'.format( ID ), fontsize = 15 )

time_text = ax.text( 0.1, 0.9, '', transform = ax.transAxes, fontsize = 13 )

ax.set_xlabel \
  ( r'$x^{{1}}\ \left[\mathrm{{{:}}}\right]$'.format( X1Units ), fontsize = 15 )
ax.set_ylabel( Field  + ' ' + r'$\left[\mathrm{{{:}}}\right]$' \
                              .format( DataUnits[2:-2] ) )

ax.set_xlim( xL, xH )
ax.set_ylim( vmin, vmax )

if UseLogScale_Y: ax.set_yscale( 'log' )
if UseLogScale_X:
    xL = max( xL, xL + 0.25 * dX10[0] )
    ax.set_xlim( xL, xH )
    ax.set_xscale( 'log' )

if PlotMesh: mesh, = ax.plot( [],[], 'b.', label = 'mesh boundaries'    )
if ShowIC: IC,     = ax.plot( [],[], 'r.', label = r'$u\left(0\right)$' )
line,              = ax.plot( [],[], 'k.', label = r'$u\left(t\right)$' )

def InitializeFrame():

    line.set_data([],[])
    time_text.set_text('')
    if ShowIC:   IC  .set_data([],[])
    if PlotMesh: mesh.set_data([],[])

    if ShowIC and PlotMesh: ret = ( line, time_text, IC, mesh )
    elif ShowIC:            ret = ( line, time_text, IC )
    elif PlotMesh:          ret = ( line, time_text, mesh )
    else:                   ret = ( line, time_text )

    return ret

def UpdateFrame( t ):

    print('    {:}/{:}'.format( t, nSS ) )
    Data, DataUnits, X1_C, dX1, Time = f(t)

    time_text.set_text( r'$t={:.3e}\ \left[\mathrm{{{:}}}\right]$' \
                        .format( Time, TimeUnits ) )

    line             .set_data( X1_C , Data .flatten() )
    if ShowIC:   IC  .set_data( X1_C0, Data0.flatten() )
    if PlotMesh: mesh.set_data( X1_C - 0.5 * dX1, \
                                0.5 * ( vmin + vmax ) \
                                  * np.ones( dX1.shape[0] ) )

    if ShowIC and PlotMesh: ret = ( line, time_text, IC, mesh )
    elif ShowIC:            ret = ( line, time_text, IC )
    elif PlotMesh:          ret = ( line, time_text, mesh )
    else:                   ret = ( line, time_text )

    return ret

ax.set_ylim( vmin, vmax )
ax.legend( prop = {'size':12} )

anim = animation.FuncAnimation( fig, UpdateFrame, \
                                init_func = InitializeFrame, \
                                frames = nSS, \
                                blit = True )

fps = max( 1, nSS / MovieRunTime )

print( '\n  Making movie' )
print( '  ------------' )
anim.save( MovieName, fps = fps, dpi = 300 )

import os
os.system( 'rm -rf __pycache__ ' )
