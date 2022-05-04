#!/usr/bin/env python3

import numpy as np
import subprocess
from os.path import isfile
import matplotlib.pyplot as plt
from matplotlib import animation

from UtilitiesModule import ChoosePlotFile, GetData, GetNorm, GetFileArray
from MakeDataFile import MakeDataFile

# --- Get user's HOME directory ---

HOME = subprocess.check_output( ["echo $HOME"], shell = True )
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# --- Get user's THORNADO_DIR directory ---

THORNADO_DIR = subprocess.check_output( ["echo $THORNADO_DIR"], shell = True )
THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

#### ========== User Input ==========

# Specify ID to be used for naming purposes
ID = 'Advection1D'

# Specify directory containing plotfiles
DataDirectory = THORNADO_DIR + 'SandBox/AMReX/Euler_Relativistic_IDEAL_MR/'

# Specify plot file base name
PlotFileBaseName = ID + '.plt'

# Specify field to plot
Field = 'PF_D'

# Specify to plot in log-scale
UseLogScale = False

# Specify whether or not to use physical units
UsePhysicalUnits = False

# Specify coordinate system (currently supports 'cartesian' and 'spherical' )
CoordinateSystem = 'cartesian'

# Specify max level of refinement (-1 is MaxLevel)
MaxLevel = -1

Verbose = True

UseCustomLimits = False
ymin = 0.0
ymax = 2.0

MovieRunTime = 10.0 # seconds

#### ====== End of User Input =======

DataFileName = '.{:s}_MovieData1D.dat'.format( ID )
MovieName    = 'mov.{:s}_{:s}.mp4'.format( ID, Field )

# Append "/" to DataDirectory, if not present
if( not DataDirectory[-1] == '/' ): DataDirectory += '/'

TimeUnit = ''
LengthUnit = ''
if UsePhysicalUnits:

    TimeUnit   = 'ms'
    LengthUnit = 'km'

Data, DataUnit, X1, X2, X3, dX1, dX2, dX3, xL, xU, nX \
  = GetData( DataDirectory, PlotFileBaseName, Field, \
             CoordinateSystem, UsePhysicalUnits, \
             MaxLevel = MaxLevel, ReturnMesh = True )

nDims = 1
if nX[1] > 1: nDims += 1
if nX[2] > 1: nDims += 1

assert ( nDims == 1 ), \
       'Invalid nDims: {:d}\nnDims must equal 1'.format( nDims )

FileArray = GetFileArray( DataDirectory, PlotFileBaseName )
#FileArray = FileArray[0:10]

MakeDataFile( Field, DataDirectory, DataFileName, \
              PlotFileBaseName, CoordinateSystem, FileArray, \
              UsePhysicalUnits = UsePhysicalUnits, \
              MaxLevel = MaxLevel )

assert isfile( DataFileName ), \
       'File: {:s} does not exist.'.format( DataFileName )

f = open( DataFileName )
header = f.readline()[16:-2]
DataUnit = f.readline()[9:-1]
DataShape = ( [ np.int64( dim ) for dim in header.split( ',' ) ] )

if UsePhysicalUnits: nOS = 12
else:                nOS = 7

Time = list( [ np.float64( t ) for t in f.readline()[nOS:-2].split(' ') ] )
Time = np.array( Time )

Data = np.loadtxt( DataFileName ).reshape( DataShape )

fig, ax = plt.subplots()

ax.set_xlim( xL[0], xU[0] )

if not UseCustomLimits:

    ymin = Data.min()
    ymax = Data.max()

ax.set_ylim( ymin, ymax )

ax.set_xlabel( 'X1' + ' ' + LengthUnit )
ax.set_ylabel( Field + ' ' + DataUnit )

Width     = xU[0] - xL[0]
Height    = ymax - ymin
time_text = plt.text( xL[0] + 0.5 * Width, ymin + 0.7 * Height, '' )

for iX1 in range( X1.shape[0] ):
    ax.axvline( X1[iX1] - 0.5 * dX1[iX1] )

if( UseLogScale ): ax.set_yscale( 'log' )

#IC,   = ax.plot([],[], color = 'red',   linewidth = 2 )
#line, = ax.plot([],[], color = 'black', linewidth = 1 )

IC,   = ax.plot([],[], 'r.' )
line, = ax.plot([],[], 'k.' )

def f( t ):
    return Data[t]

def InitializeFrame():

    IC.set_data([],[])
    line.set_data([],[])
    time_text.set_text('')

    ret = ( line, time_text, IC )
    return ret

def UpdateFrame( t ):

    IC.set_data( X1, f(0) )
    y = Data[t]
    line.set_data( X1, y )
    time_text.set_text( 'time = {:.3e} {:}'.format( Time[t], TimeUnit ) )

    ret = ( line, time_text, IC )
    return ret

nFrames = FileArray.shape[0]

fps = max( 1, nFrames / MovieRunTime )

anim = animation.FuncAnimation( fig, UpdateFrame, \
                                init_func = InitializeFrame, \
                                frames = nFrames, \
                                blit = True )

anim.save( MovieName, fps = fps )

import os
os.system( 'rm -rf __pycache__ ' )
