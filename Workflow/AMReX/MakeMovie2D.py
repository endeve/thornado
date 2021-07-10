#!/usr/bin/env python3

import numpy as np
import subprocess
from os.path import isfile
import matplotlib.pyplot as plt
from matplotlib import animation

from UtilitiesModule import ChoosePlotFile, GetData, GetNorm
from MakeDataFile import MakeDataFile

# --- Get user's HOME directory ---

HOME = subprocess.check_output( ["echo $HOME"], shell = True )
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# --- Get user's THORNADO_DIR directory ---

THORNADO_DIR = subprocess.check_output( ["echo $THORNADO_DIR"], shell = True )
THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

#### ========== User Input ==========

# Specify ID to be used for naming purposes
ID = 'SAS2D'

# Specify directory containing plotfiles
DataDirectory = THORNADO_DIR + 'SandBox/AMReX/'

# Specify plot file base name
PlotFileBaseName = 'plt'

# Specify field to plot
Field = 'PF_D'

# Specify to plot in log-scale
UseLogScale = True

# Specify whether or not to use physical units
UsePhysicalUnits = True

# Specify coordinate system (currently supports 'cartesian' and 'spherical' )
CoordinateSystem = 'spherical'

# Specify colormap
cmap = 'cividis'

Verbose = True

UseCustomLimits = False
vmin = 0.0
vmax = 2.0

MovieRunTime = 10.0 # seconds

zAxisVertical = False

#### ====== End of User Input =======

DataFileName = '.{:s}_{:s}_MovieData2D.dat'.format( ID, Field )
MovieName    = 'mov.{:s}_{:s}.mp4'.format( ID, Field )

# Append "/" to DataDirectory, if not present
if( not DataDirectory[-1] == '/' ): DataDirectory += '/'

File, FileArray \
  = ChoosePlotFile( DataDirectory, PlotFileBaseName, \
                    ReturnFileArray = True, Verbose = Verbose )

TimeUnit = ''
LengthUnit = ''
if UsePhysicalUnits:

    TimeUnit   = 'ms'
    LengthUnit = 'km'

Data, DataUnit, X1, X2, X3, xL, xU, nX, Time \
  = GetData( DataDirectory, PlotFileBaseName, Field, \
             CoordinateSystem, UsePhysicalUnits, \
             ReturnTime = True, ReturnMesh = True, \
             Verbose = False )

nDims = 1
if nX[1] > 1: nDims += 1
if nX[2] > 1: nDims += 1

assert ( nDims == 2 ), \
       'Invalid nDims: {:d}\nnDims must equal 2'.format( nDims )

## To make lineout plot
## From: https://yt-project.org/doc/visualizing/
##       manual_plotting.html#line-plots
#
#oray = ds.ortho_ray( axis = 0, coords = (0,0) )
#
#plt.plot( X1, oray[Field], 'k-' )
#
#if( UseLogScale ): plt.yscale( 'log' )
#plt.show()
#plt.close()
#exit()

MakeDataFile( Field, DataDirectory, DataFileName, \
              PlotFileBaseName, CoordinateSystem, \
              UsePhysicalUnits = UsePhysicalUnits, \
              WriteExtras = False, Verbose = False )

if( not isfile( DataFileName ) ):

    print( 'File: {:s} does not exist.'.format( DataFileName ) )
    exit( 'Exiting...' )

assert isfile( DataFileName ), \
       'File: {:s} does not exist.'.format( DataFileName )

f = open( DataFileName  )
header = f.readline()[16:-2]
DataUnit = f.readline()[9:-1]
DataShape = ( [ np.int64( dim ) for dim in header.split( ',' ) ] )

if UsePhysicalUnits:

    Time = list( [ np.float64( t ) \
                 for t in f.readline()[12:-2].split(' ') ] )

else:

    Time = list( [ np.float64( t ) \
                 for t in f.readline()[7:-2].split(' ') ] )

Time = np.array( Time )

Data = np.loadtxt( DataFileName ).reshape( DataShape )

if CoordinateSystem == 'spherical':

    polar = True

    figsize = (16,9)

else:

    polar = False

    figsize = (8,6)

fig = plt.figure( figsize = figsize )
ax  = fig.add_subplot( 111, polar = polar )

def f(t):
    return Data[t]

if not UseCustomLimits:

    vmin = Data.min()
    vmax = Data.max()

Norm = GetNorm( UseLogScale, Data, vmin = vmin, vmax = vmax )

if CoordinateSystem == 'spherical':

    im = ax.pcolormesh( X2, X1, f(0)[:,:], \
                        cmap = cmap, \
                        norm = Norm, \
                        shading = 'nearest' )

    ax.set_thetamin( 0.0  )
    ax.set_thetamax( 180.0 )
    ax.set_theta_direction( -1 )
    ax.set_theta_zero_location( 'W' ) # z-axis horizontal

    if( zAxisVertical ):
        ax.set_theta_zero_location( 'N' ) # z-axis vertical
        time_text = plt.text( 0.5 * np.pi / 2, xU[0] * ( 1.0 + 0.3 ), '' )
        cbar = fig.colorbar( im )
    else:
        ax.set_theta_zero_location( 'W' ) # z-axis horizontal
        time_text = plt.text( 0.9 *np.pi / 2, xU[0] * ( 1.0 + 0.3 ), '' )
        ax.set_position( [0.1,-0.45,0.7,2] )
        cax = fig.add_axes( [0.85,0.1,0.03,0.8] )
        cbar = fig.colorbar( im, cax = cax )

else:

    im = ax.pcolormesh( X1, X2, f(0)[:,:], \
                        cmap = cmap, \
                        norm = Norm, \
                        shading = 'nearest' )

    ax.set_xlim( xL[0], xU[0] )
    ax.set_ylim( xL[1], xU[1] )

    Width  = xU[0] - xL[0]
    Height = xU[1] - xL[1]

    time_text \
      = plt.text( xL[0] + 0.5 * Width, xL[1] + 0.9 * Height, '' )

    cbar = fig.colorbar( im )

cbar.set_label( Field + ' ' + DataUnit )

def UpdateFrame(t):

    im.set_array( f(t)[:,:].flatten() )
    time_text.set_text( 'Time = {:.3e} {:}'.format( Time[t], TimeUnit ) )

    ret = ( im, time_text )
    return ret

nFrames = FileArray.shape[0]

# Call the animator
anim = animation.FuncAnimation \
         ( fig, UpdateFrame, frames = nFrames, \
           blit = True)

fps = max( 1, nFrames / MovieRunTime )

anim.save( MovieName, fps = fps )

import os
os.system( 'rm -rf __pycache__ ' )
