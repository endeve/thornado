#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

from UtilitiesModule import GetData, GetNorm, GetFileArray
from MakeDataFile import MakeDataFile, ReadHeader

#### ========== User Input ==========

# ID to be used for naming purposes
ID = 'KelvinHelmholtz2D'

# Directory containing AMReX plotfiles
PlotFileDirectory \
  = 'thornado/SandBox/AMReX/Euler_Relativistic_IDEAL_MR/'

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

# Colormap
cmap = 'viridis'

# First and last snapshots and number of snapshots to include in movie
SSi = -1 # -1 -> SSi = 0
SSf = -1 # -1 -> PlotFileArray.shape[0] - 1
nSS = -1 # -1 -> PlotFileArray.shape[0]

# Max level of refinement to include
MaxLevel = -1 # -1 -> use all levels

Verbose = True

UseCustomLimits = False
vmin = 0.0
vmax = 2.0

MovieRunTime = 10.0 # seconds

zAxisVertical = False

#### ====== End of User Input =======

DataFileDirectory = '.{:s}_{:s}_MovieData2D'.format( ID, Field )
MovieName         = 'mov.{:s}_{:s}.mp4'.format( ID, Field )

# Append "/" if not present
if( not PlotFileDirectory[-1] == '/' ): PlotFileDirectory += '/'
if( not DataFileDirectory[-1] == '/' ): DataFileDirectory += '/'

TimeUnit = ''
if UsePhysicalUnits: TimeUnit = 'ms'

MakeDataFile( Field, PlotFileDirectory, DataFileDirectory, \
              PlotFileBaseName, CoordinateSystem, \
              SSi = SSi, SSf = SSf, nSS = nSS, \
              UsePhysicalUnits = UsePhysicalUnits, \
              MaxLevel = MaxLevel, Verbose = Verbose )

PlotFileArray \
  = GetFileArray( PlotFileDirectory, PlotFileBaseName, \
                  SSi = SSi, SSf = SSf, nSS = nSS )

def f(t):

    DataFile = DataFileDirectory + PlotFileArray[t] + '.dat'

    DataShape, DataUnits, Time, X1_C, X2_C, X3_C, dX1, dX2, dX3 \
      = ReadHeader( DataFile )

    Data = np.loadtxt( DataFile ).reshape( np.int64( DataShape ) )

    return Data, DataUnits, X1_C, X2_C, dX1, dX2, Time

Data, DataUnits, \
  X1_C, X2_C, X3_C, dX1, dX2, dX3, xL, xU, nX \
    = GetData( PlotFileDirectory, PlotFileBaseName, Field, \
               CoordinateSystem, UsePhysicalUnits, \
               argv = [ 'a' ], \
               MaxLevel = MaxLevel, \
               ReturnTime = False, ReturnMesh = True )

xL = np.array( [ X1_C[0 ]-0.5*dX1[0 ], X2_C[0 ]-0.5*dX2[0 ] ], np.float64 )
xU = np.array( [ X1_C[-1]+0.5*dX1[-1], X2_C[-1]+0.5*dX2[-1] ], np.float64 )

if not UseCustomLimits:
    vmin = Data.min()
    vmax = Data.max()

Norm = GetNorm( UseLogScale, Data, vmin = vmin, vmax = vmax )

X1v, X2v = np.meshgrid( X1_C, X2_C, indexing = 'ij' )

fig = plt.figure()
ax  = fig.add_subplot( 111 )

im = ax.pcolormesh( X1v, X2v, Data, \
                    cmap = cmap, \
                    norm = Norm, \
                    shading = 'nearest' )

ax.set_xlim( xL[0], xU[0] )
ax.set_ylim( xL[1], xU[1] )
#ax.set_aspect( 'auto' )

time_text = plt.text( 0.5, 0.9, '', transform = ax.transAxes )

cbar = fig.colorbar( im )

cbar.set_label( Field + ' ' + DataUnits )

def UpdateFrame(t):

    Data, DataUnits, X1_C, X2_C, dX1, dX2, Time = f(t)

    X1v, X2v = np.meshgrid( X1_C, X2_C, indexing = 'ij' )

    im = ax.pcolormesh( X1v, X2v, Data, \
                        cmap = cmap, \
                        norm = Norm, \
                        shading = 'nearest' )

    im.set_array( Data.flatten() )
    time_text.set_text( 'Time = {:.3e} {:}'.format( Time, TimeUnit ) )

    ret = ( im, time_text )
    return ret

if nSS < 0: nSS = PlotFileArray.shape[0]

# Call the animator
anim = animation.FuncAnimation \
         ( fig, UpdateFrame, frames = nSS, \
           blit = True)

fps = max( 1, nSS / MovieRunTime )

anim.save( MovieName, fps = fps )

import os
os.system( 'rm -rf __pycache__ ' )
