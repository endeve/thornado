#!/usr/bin/env python3

import numpy as np
import subprocess
from os.path import isfile
import matplotlib.pyplot as plt
from matplotlib import animation

from UtilitiesModule import GetNorm, GetFileArray
from MakeDataFile import MakeDataFile, ReadHeader

# --- Get user's HOME directory ---

HOME = subprocess.check_output( ["echo $HOME"], shell = True )
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# --- Get user's THORNADO_DIR directory ---

THORNADO_DIR = subprocess.check_output( ["echo $THORNADO_DIR"], shell = True )
THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

#### ========== User Input ==========

# ID to be used for naming purposes
ID = 'RiemannProblem1D'

# Directory containing AMReX plotfiles
PlotFileDirectory \
  = THORNADO_DIR + 'SandBox/AMReX/Euler_Relativistic_IDEAL_MR/'

# PlotFile base name (e.g., Advection2D.plt######## -> Advection2D.plt )
PlotFileBaseName = ID + '.plt'

# Field to plot
Field = 'PF_D'

# Plot data in log10-scale?
UseLogScale = False

# Unit system of the data
UsePhysicalUnits = False

# Coordinate system (currently supports 'cartesian' and 'spherical' )
CoordinateSystem = 'cartesian'

# First and last snapshots and number of snapshots to include in movie
SSi = -1 # -1 -> SSi = 0
SSf = -1 # -1 -> PlotFileArray.shape[0] - 1
nSS = -1 # -1 -> PlotFileArray.shape[0]

# Max level of refinement to include
MaxLevel = -1

# Include initial conditions in movie?
ShowIC = True

PlotMesh = False

Verbose = True

UseCustomLimits = False
ymin = 0.0
ymax = 2.0

MovieRunTime = 10.0 # seconds

#### ====== End of User Input =======

DataFileDirectory = '.{:s}_{:s}_MovieData1D'.format( ID, Field )
MovieName         = 'mov.{:s}_{:s}.mp4'.format( ID, Field )

# Append "/" if not present
if( not PlotFileDirectory[-1] == '/' ): PlotFileDirectory += '/'
if( not DataFileDirectory[-1] == '/' ): DataFileDirectory += '/'

MakeDataFile( Field, PlotFileDirectory, DataFileDirectory, \
              PlotFileBaseName, CoordinateSystem, \
              SSi = SSi, SSf = SSf, nSS = nSS, \
              UsePhysicalUnits = UsePhysicalUnits, \
              MaxLevel = MaxLevel, Verbose = Verbose )

PlotFileArray = GetFileArray( PlotFileDirectory, PlotFileBaseName )

if SSi < 0: SSi = 0
if SSf < 0: SSf = PlotFileArray.shape[0] - 1
if nSS < 0: nSS = PlotFileArray.shape[0]

fig = plt.figure()
ax  = fig.add_subplot( 111 )

if UseCustomLimits: ax.set_ylim( ymin, ymax )
ymin = +np.inf
ymax = -np.inf
for j in range( nSS ):
    iSS = SSi + np.int64( ( SSf - SSi ) / ( nSS - 1 ) * j )
    DataFile = DataFileDirectory + PlotFileArray[iSS] + '.dat'
    Data = np.loadtxt( DataFile )
    ymin = min( ymin, Data.min() )
    ymax = max( ymax, Data.max() )
if not UseCustomLimits: ax.set_ylim( ymin, ymax )

def f(t):

    iSS = SSi + np.int64( ( SSf - SSi + 1 ) / nSS ) * t

    DataFile = DataFileDirectory + PlotFileArray[iSS] + '.dat'

    DataShape, DataUnits, Time, X1_C, X2_C, X3_C, dX1, dX2, dX3 \
      = ReadHeader( DataFile )

    Data = np.loadtxt( DataFile ).reshape( np.int64( DataShape ) )

    global ymin, ymax
    ymin = min( ymin, Data.min() )
    ymax = max( ymax, Data.max() )

    return Data, DataUnits, X1_C, dX1, Time

Data0, DataUnits, X1_C0, dX10, Time0 = f(0)

xL = np.array( [ X1_C0[0 ]-0.5*dX10[0 ] ], np.float64 )
xU = np.array( [ X1_C0[-1]+0.5*dX10[-1] ], np.float64 )

Norm = GetNorm( UseLogScale, Data0, vmin = ymin, vmax = ymax )

TimeUnits = ''
LengthUnits = ''
if UsePhysicalUnits:
    TimeUnits   = 'ms'
    LengthUnits = 'km'

ax.set_xlim( xL[0], xU[0] )

ax.set_xlabel( 'X1' + ' ' + LengthUnits )
ax.set_ylabel( Field + ' ' + DataUnits )

time_text = plt.text( 0.5, 0.7, '', transform = ax.transAxes )

if( UseLogScale ): ax.set_yscale( 'log' )

line, = ax.plot( [],[], 'k.', label = 'u(t)' )
if ShowIC: IC, = ax.plot( [],[], 'r.', label = 'u(0)' )
if PlotMesh: mesh, = ax.plot( [],[], 'b.', label = 'mesh' )

def InitializeFrame():

    line.set_data([],[])
    time_text.set_text('')
    if ShowIC: IC.set_data([],[])
    if PlotMesh: mesh.set_data([],[])

    if ShowIC and PlotMesh: ret = ( line, time_text, IC, mesh )
    elif ShowIC:            ret = ( line, time_text, IC )
    elif PlotMesh:          ret = ( line, time_text, mesh )
    else:                   ret = ( line, time_text )

    return ret

def UpdateFrame( t ):

    Data, DataUnits, X1_C, dX1, Time = f(t)

    line.set_data( X1_C, Data.flatten() )
    time_text.set_text( 'Time = {:.3e} {:}'.format( Time, TimeUnits ) )
    if ShowIC: IC.set_data( X1_C0, Data0.flatten() )
    if PlotMesh: mesh.set_data( X1_C - 0.5 * dX1, 0.5 * ( ymin + ymax ) )

    if ShowIC and PlotMesh: ret = ( line, time_text, IC, mesh )
    elif ShowIC:            ret = ( line, time_text, IC )
    elif PlotMesh:          ret = ( line, time_text, mesh )
    else:                   ret = ( line, time_text )

    return ret

ax.set_ylim( ymin, ymax )
ax.legend()

anim = animation.FuncAnimation( fig, UpdateFrame, \
                                init_func = InitializeFrame, \
                                frames = nSS, \
                                blit = True )


fps = max( 1, nSS / MovieRunTime )

anim.save( MovieName, fps = fps )

import os
os.system( 'rm -rf __pycache__ ' )
