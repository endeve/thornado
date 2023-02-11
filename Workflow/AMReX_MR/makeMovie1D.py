#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

from UtilitiesModule import GetNorm
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

# Directory containing AMReX Plotfiles
PlotfileDirectory = 'thornado/SandBox/AMReX/'
#PlotfileDirectory \
#  = '/home/kkadoogan/Work/Codes/thornado/\
#SandBox/AMReX/Euler_Relativistic_IDEAL_MR/'

# Plotfile base name (e.g., Advection1D.plt######## -> Advection1D.plt )
PlotfileBaseName = ID + '.plt'

# Field to plot
Field = 'PF_D'

# Plot data in log10-scale?
UseLogScale_Y = False
UseLogScale_X = False

# Unit system of the data
UsePhysicalUnits = False

# Coordinate system (currently supports 'cartesian' and 'spherical' )
CoordinateSystem = 'cartesian'

# First and last snapshots and number of snapshots to include in movie
SSi = -1 # -1 -> SSi = 0
SSf = -1 # -1 -> PlotfileArray.shape[0] - 1
nSS = -1 # -1 -> PlotfileArray.shape[0]

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
if( not PlotfileDirectory[-1] == '/' ): PlotfileDirectory += '/'
if( not DataDirectory    [-1] == '/' ): DataDirectory     += '/'

PlotfileArray \
  = MakeDataFile( Field, PlotfileDirectory, DataDirectory, \
                  PlotfileBaseName, CoordinateSystem, \
                  SSi = SSi, SSf = SSf, nSS = nSS, \
                  UsePhysicalUnits = UsePhysicalUnits, \
                  MaxLevel = MaxLevel, forceChoice = False, ow = False, \
                  Verbose = Verbose )

if nSS < 0: nSS = PlotfileArray.shape[0]

fig = plt.figure()
ax  = fig.add_subplot( 111 )

if not UseCustomLimits:
    vmin = +np.inf
    vmax = -np.inf
    for j in range( nSS ):
        DataFile \
          = DataDirectory + PlotfileArray[j] + '/{:}.dat'.format( Field )
        DataShape, DataUnits, MinVal, MaxVal = ReadHeader( DataFile )
        vmin = min( vmin, MinVal )
        vmax = max( vmax, MaxVal )
ax.set_ylim( vmin, vmax )

def f(t):

    FileDirectory = DataDirectory + PlotfileArray[t] + '/'

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

Data0, DataUnits, X1_C0, dX1, Time = f(0)

xL = X1_C0[0 ] - 0.5 * dX1[0 ]
xH = X1_C0[-1] + 0.5 * dX1[-1]

Norm = GetNorm( UseLogScale_Y, Data0, vmin = vmin, vmax = vmax )

TimeUnits = ''
LengthUnits = '[]'
if UsePhysicalUnits:
    TimeUnits   = 'ms'
    LengthUnits = '[km]'

ax.set_xlim( xL, xH )

ax.set_xlabel( 'X1_C' + ' ' + LengthUnits )
ax.set_ylabel( Field  + ' ' + DataUnits   )

time_text = plt.text( 0.5, 0.7, '', transform = ax.transAxes )

if( UseLogScale_Y ): ax.set_yscale( 'log' )
if( UseLogScale_X ): ax.set_xscale( 'log' )

if PlotMesh: mesh, = ax.plot( [],[], 'b.', label = 'mesh boundaries' )
if ShowIC: IC,     = ax.plot( [],[], 'r.', label = 'u(0)' )
line,              = ax.plot( [],[], 'k.', label = 'u(t)' )

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

    Data, DataUnits, X1_C, dX1, Time = f(t)

    time_text.set_text( 'Time = {:.3e} {:}'.format( Time, TimeUnits ) )

    line             .set_data( X1_C , Data .flatten() )
    if ShowIC:   IC  .set_data( X1_C0, Data0.flatten() )
    if PlotMesh: mesh.set_data( X1_C - 0.5 * dX1, 0.5 * ( vmin + vmax ) )

    if ShowIC and PlotMesh: ret = ( line, time_text, IC, mesh )
    elif ShowIC:            ret = ( line, time_text, IC )
    elif PlotMesh:          ret = ( line, time_text, mesh )
    else:                   ret = ( line, time_text )

    return ret

ax.set_ylim( vmin, vmax )
ax.legend()

anim = animation.FuncAnimation( fig, UpdateFrame, \
                                init_func = InitializeFrame, \
                                frames = nSS, \
                                blit = True )

fps = max( 1, nSS / MovieRunTime )

anim.save( MovieName, fps = fps )

import os
os.system( 'rm -rf __pycache__ ' )
