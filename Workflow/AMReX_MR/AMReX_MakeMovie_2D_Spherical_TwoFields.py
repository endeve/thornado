#!/usr/bin/env python3

"""
Generate a movie showing two fields side-by-side
from data files created in MakeDataFiles.py.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import subprocess

# --- Get user's HOME directory ---
HOME = subprocess.check_output( ["echo $HOME"], shell = True)
HOME = HOME[:-1].decode( "utf-8" ) + '/'

############# User input #############

# First

Root1   = 'M1.4_Rs180_Mdot0.3'
suffix1 = ''

DataDirectory1 = HOME + '{:}/'.format( Root1 )

PlotFileBaseName1 = 'plt_{:}{:}'.format( Root1, suffix1 )

Field1 = 'PolytropicConstant' # Field to plot

cmap1 = 'Purples' # Color scheme

UseLogScale1 = True # Use log scale?

# Second

Root2   = 'M1.4_Rs180_Mdot0.3'
suffix2 = ''

DataDirectory2 = HOME + '{:}/'.format( Root2 )

PlotFileBaseName2 = 'plt_{:}{:}'.format( Root2, suffix2 )

Field2 = 'PolytropicConstant' # Field to plot

cmap2 = 'Purples' # Color scheme

UseLogScale2 = True # Use log scale?

# Global

MovieName = 'SASI_{:}_{:}'.format( Field1, Field2 )

MovieRunTime = 30.0 # seconds

UsePhysicalUnits = True # Does your data use physical units?

UseCustomTicks = True # Define limits for colorbar near line 70

zAxisVertical = True # Orient z-axis

############# End of user input #############

ID1           = '{:}{:}_{:}'.format( Root1, suffix1, Field1 )
DataFileName1 = '{:}_MovieData.dat'.format( ID1 )
TimeFileName1 = '{:}_MovieTime.dat'.format( ID1 )

ID2           = '{:}{:}_{:}'.format( Root2, suffix2, Field2 )
DataFileName2 = '{:}_MovieData.dat'.format( ID2 )
TimeFileName2 = '{:}_MovieTime.dat'.format( ID2 )

from MakeDataFiles import *

xL, xH, nX, FileArray1 \
  = MakeDataFile( Field1, DataDirectory1, DataFileName1, TimeFileName1, \
                  PlotFileBaseName1, UsePhysicalUnits )

xL, xH, nX, FileArray2 \
  = MakeDataFile( Field2, DataDirectory2, DataFileName2, TimeFileName2, \
                  PlotFileBaseName2, UsePhysicalUnits )

print( 'Reading in data files...' )

Data1 = np.loadtxt( DataFileName1 ).reshape( \
                                    (FileArray1.shape[0],nX[0],nX[1]) )
Data2 = np.loadtxt( DataFileName2 ).reshape( \
                                    (FileArray2.shape[0],nX[0],nX[1]) )

Time = np.loadtxt( TimeFileName1 ) # These should be the same for both fields

fig = plt.figure( figsize = (8,6) )
ax  = fig.add_subplot( 111, polar = True )
xL  = xL.to_ndarray()
xH  = xH.to_ndarray()
X1  = np.linspace( xL[0], xH[0], nX[0] )
X2  = np.linspace( xL[1], xH[1], nX[1] )
theta1, r = np.meshgrid( X2, X1 )
theta2, r = np.meshgrid( 2.0 * np.pi - X2, X1 )

if( UseCustomTicks ):

    vmin1 = min( +np.inf, np.min( Data1 ) )
    vmax1 = max( -np.inf, np.max( Data1 ) )
    vmin2 = min( +np.inf, np.min( Data2 ) )
    vmax2 = max( -np.inf, np.max( Data2 ) )

    if( UseLogScale1 ):

        if  ( vmax1 < 0.0 ):
            ticks1 = np.linspace( -np.log10(-vmin1), -np.log10(-vmax1), 5 )
        elif( vmin1 < 0.0 ):
            ticks1 = np.linspace( -np.log10(-vmin1), +np.log10(+vmax1), 5 )
        else:
            ticks1 = np.logspace( +np.log10(+vmin1), +np.log10(+vmax1), 5 )

    else:

        ticks1 = np.linspace( vmin1, vmax1, 5 )

    if( UseLogScale2 ):

        if  ( vmax2 < 0.0 ):
            ticks2 = np.linspace( -np.log10(-vmin2), -np.log10(-vmax2), 5 )
        elif( vmin2 < 0.0 ):
            ticks2 = np.linspace( -np.log10(-vmin2), +np.log10(+vmax2), 5 )
        else:
            ticks2 = np.logspace( +np.log10(+vmin2), +np.log10(+vmax2), 5 )

    else:

        ticks2 = np.linspace( vmin2, vmax2, 5 )

    ticklabels1 = []
    ticklabels2 = []
    for tick in ticks1:
        ticklabels1.append( '{:.3e}'.format( tick ) )
    for tick in ticks2:
        ticklabels2.append( '{:.3e}'.format( tick ) )

else:

    vmin1 = min( +np.inf, np.min( Data1 ) )
    vmax1 = max( -np.inf, np.max( Data1 ) )
    vmin2 = min( +np.inf, np.min( Data2 ) )
    vmax2 = max( -np.inf, np.max( Data2 ) )

if( UseLogScale1 ):

    from matplotlib.colors import LogNorm, SymLogNorm

    if( np.any( Data1 < 0.0 ) ):
        Norm1 = SymLogNorm( vmin = vmin1, vmax = vmax1, \
                            linthresh = 1.0e2, base = 10 )
    else:
        Norm1 = LogNorm   ( vmin = vmin1, vmax = vmax1 )

else:

    Norm1 = plt.Normalize ( vmin = vmin1, vmax = vmax1 )

if( UseLogScale2 ):

    from matplotlib.colors import LogNorm, SymLogNorm

    if( np.any( Data2 < 0.0 ) ):
        Norm2 = SymLogNorm( vmin = vmin2, vmax = vmax2, \
                            linthresh = 1.0e2, base = 10 )
    else:
        Norm2 = LogNorm   ( vmin = vmin2, vmax = vmax2 )

else:

    Norm2 = plt.Normalize ( vmin = vmin2, vmax = vmax2 )

def f1(t):
    return Data1[t]
def f2(t):
    return Data2[t]

# Taken from:
# https://brushingupscience.com/2016/06/21/matplotlib-animations-the-easy-way/
im1 = ax.pcolormesh( theta1, r, f1(0)[:-1,:-1], \
                     cmap = cmap1, \
                     norm = Norm1 )
im2 = ax.pcolormesh( theta2, r, f2(0)[:-1,:-1], \
                     cmap = cmap2, \
                     norm = Norm2 )

# Limits on coordinate axes

ax.set_thetamin( 0.0   )
ax.set_thetamax( 360.0 )
ax.set_theta_direction( -1 )
ax.set_xticklabels( [] )

if( zAxisVertical ):
    ax.set_theta_zero_location( 'N' ) # z-axis vertical
    time_text = plt.text( 0.1 * np.pi / 2.0, xH[0] * ( 1.0 + 0.1 ), '' )
else:
    ax.set_theta_zero_location( 'W' ) # z-axis horizontal
    time_text = plt.text( 0.9 * np.pi / 2.0, xH[0] * ( 1.0 + 0.3 ), '' )

cb1axes = fig.add_axes( [ 0.81,  0.1, 0.03, 0.8 ] )
cb2axes = fig.add_axes( [ 0.185, 0.1, 0.03, 0.8 ] )
if( UseCustomTicks ):
    cbar1 = fig.colorbar( im1, cax = cb1axes, ticks = ticks1 )
    cbar1.ax.set_yticklabels( ticklabels1 )
    cbar2 = fig.colorbar( im2, cax = cb2axes, ticks = ticks2 )
    cbar2.ax.set_yticklabels( ticklabels2 )
else:
    cbar1 = fig.colorbar( im1, cax = cb1axes )
    cbar2 = fig.colorbar( im2, cax = cb2axes )
cbar1.ax.set_ylabel( Field1 )
cbar2.ax.set_ylabel( Field2 )

cb2axes.yaxis.set_ticks_position( 'left' )
cb2axes.yaxis.set_label_position( 'left' )

TimeUnit = ''
if( UsePhysicalUnits ): TimeUnit = ' ms'

def UpdateFrame(t):
    im1.set_array( f1(t)[:-1,:-1].flatten() )
    im2.set_array( f2(t)[:-1,:-1].flatten() )
    time_text.set_text('time = {:d}{:}'.format( np.int( Time[t] ), TimeUnit ) )
    return im1, im2,

# Call the animator

print( 'Making movie...' )

nFrames = max( FileArray1.shape[0], FileArray2.shape[0] )
fps = nFrames / MovieRunTime

anim \
  = animation.FuncAnimation \
      ( fig, UpdateFrame, frames = nFrames, blit = True )

anim.save( '{:}_Movie.mp4'.format( MovieName ), fps = fps )
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
