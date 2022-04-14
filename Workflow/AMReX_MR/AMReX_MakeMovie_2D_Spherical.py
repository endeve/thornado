#!/usr/bin/env python3

"""
Generate a movie from data files created in MakeDataFiles.py.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import subprocess

# --- Get user's HOME directory ---
HOME = subprocess.check_output( ["echo $HOME"], shell = True)
HOME = HOME[:-1].decode( "utf-8" ) + '/'

############# User input #############

Root   = 'M1.4_Rs180_Mdot0.3'
suffix = ''

DataDirectory = HOME + '{:}/'.format( Root )

PlotFileBaseName = 'plt_{:}{:}'.format( Root, suffix )

Field = 'PF_D' # Field to plot

cmap = 'Purples' # Color scheme for movie

UseLogScale = True # Use log scale for field?

MovieName = 'SASI_{:}'.format( Field )

MovieRunTime = 30.0 # seconds

UsePhysicalUnits = True # Does your data use physical units?

UseCustomTicks = True # Define limits for colorbar near line 70

zAxisVertical = False # Orient z-axis

WriteExtras = False # Set to true if generating data on
                    # external machine (Summit, ACF, etc.)

############# End of user input #############

ID = '{:}{:}_{:}'.format( Root, suffix, Field )

DataFileName = '{:}_MovieData.dat'.format( ID )
TimeFileName = '{:}_MovieTime.dat'.format( ID )

from MakeDataFiles import *

xL, xH, nX, FileArray \
  = MakeDataFile( Field, DataDirectory, DataFileName, TimeFileName, \
                  PlotFileBaseName, UsePhysicalUnits, \
                  WriteExtras )

print( 'Reading in data file...' )
Data = np.loadtxt( DataFileName ).reshape( \
                                    (FileArray.shape[0],nX[0],nX[1]) )
Time = np.loadtxt( TimeFileName )

fig = plt.figure( figsize = (8,6) )
ax  = fig.add_subplot( 111, polar = True )
xL  = xL.to_ndarray()
xH  = xH.to_ndarray()
X1  = np.linspace( xL[0], xH[0], nX[0]+1 )
X2  = np.linspace( xL[1], xH[1], nX[1]+1 )
theta, r = np.meshgrid( X2, X1 )

vmin = min( +np.inf, np.min( Data ) )
vmax = max( -np.inf, np.max( Data ) )

if( UseCustomTicks ):

    vmin = min( +np.inf, Data.min() )
    vmax = max( -np.inf, Data.max() )

    if( UseLogScale ):

        if  ( vmax < 0.0 ):
            ticks = np.linspace( -np.log10(-vmin), -np.log10(-vmax), 5 )
        elif( vmin < 0.0 ):
            ticks = np.linspace( -np.log10(-vmin), +np.log10(+vmax), 5 )
        else:
            ticks = np.logspace( +np.log10(+vmin), +np.log10(+vmax), 5 )

    else:

        ticks = np.linspace( vmin, vmax, 5 )

    ticklabels = []
    for tick in ticks:
        ticklabels.append( '{:.3e}'.format( tick ) )

if( UseLogScale ):

    from matplotlib.colors import LogNorm, SymLogNorm

    if( np.any( Data < 0.0 ) ):

        Norm = SymLogNorm( vmin = vmin, vmax = vmax, \
                           linthresh = 1.0e2, base = 10 )

    else:

        Norm = LogNorm( vmin = vmin, vmax = vmax )
else:

    Norm = plt.Normalize( vmin = vmin, vmax = vmax )

def f(t):
    return Data[t]

# Taken from:
# https://brushingupscience.com/2016/06/21/matplotlib-animations-the-easy-way/
im = ax.pcolormesh( theta, r, f(0)[:,:], \
                    cmap = cmap, \
                    norm = Norm )

# Limits on coordinate axes

ax.set_thetamin( 180.0/np.pi * X2[0]  )
ax.set_thetamax( 180.0/np.pi * X2[-1] )
ax.set_theta_direction( -1 )

if( zAxisVertical ):
    ax.set_theta_zero_location( 'N' ) # z-axis vertical
    time_text = plt.text( 0.5 * np.pi / 2, xH[0] * ( 1.0 + 0.3 ), '' )
else:
    ax.set_theta_zero_location( 'W' ) # z-axis horizontal
    time_text = plt.text( 0.9 *np.pi / 2, xH[0] * ( 1.0 + 0.3 ), '' )

if( UseCustomTicks ):
    cbar = fig.colorbar( im, ticks = ticks )
    cbar.ax.set_yticklabels( ticklabels )
else:
    cbar = fig.colorbar( im )

TimeUnit = ''
if( UsePhysicalUnits ): TimeUnit = ' ms'

def UpdateFrame(t):
    im.set_array( f(t)[:,:].flatten() )
    time_text.set_text('time = {:d}{:}'.format( np.int( Time[t] ), TimeUnit ) )
    return im,

# Call the animator

print( 'Making movie...' )

nFrames = FileArray.shape[0]
fps = nFrames / MovieRunTime

anim \
  = animation.FuncAnimation \
      ( fig, UpdateFrame, frames = nFrames, blit = True )

anim.save( '{:}_Movie.mp4'.format( MovieName ), fps = fps )
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
