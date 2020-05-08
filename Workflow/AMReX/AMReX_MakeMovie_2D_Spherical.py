#!/usr/local/bin/python3

"""
Generate a movie from data files created in MakeDataFiles.py.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

############# User input #############
ProblemName = 'SASI' # Only used for name of movie file

# Define custom variables as seen near line 40
Field = 'PF_D'

DataFileName = 'MovieData_{:}.dat'.format( Field )

UseLogScale = True # Do you want your field plotted in log scale?

cmap = 'Purples' # Color scheme for movie

UseCustomTicks = True # Define limits for colorbar near line 166

zAxisVertical = False # Orient z-axis
############# End of user input #############

from MakeDataFiles import *

MakeDataFile( Field, DataFileName )

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

if( UseLogScale ):

    if( np.any( Data < 0.0 ) ):
        def f(t):
            return np.sign( Data[t] ) * np.log10( np.abs( Data[t] ) + 1.0 )
    else:
        def f(t):
            return np.log10( Data[t] )

else:

    def f(t):
        return Data[t]

if( UseCustomTicks ):

    vmin = min( +np.inf, np.min( Data ) )
    vmax = max( -np.inf, np.max( Data ) )

    if( UseLogScale ):

        vmax = np.log10( vmax )

        if( np.any( Data < 0.0 ) ):
            vmin = -np.log10( abs( np.min( Data ) ) )
        else:
          vmin = np.log10( vmin )

        ticks = np.linspace( vmin, vmax, 5 )

    else:
        ticks = np.linspace( vmin, vmax, 5 )

    ticklabels = []
    for tick in ticks:
        ticklabels.append( '{:.3e}'.format( tick ) )

else:

    vmin = min( +np.inf, np.min( Data ) )
    vmax = max( -np.inf, np.max( Data ) )

# Taken from:
# https://brushingupscience.com/2016/06/21/matplotlib-animations-the-easy-way/
im = ax.pcolormesh( theta, r, f(0)[:,:], \
                    cmap = cmap, \
                    vmin = vmin, vmax = vmax, \
                    norm = None )
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

def UpdateFrame(t):
    im.set_array( f(t)[:,:].flatten() )
    if( UsePhysicalUnits ):
      time_text.set_text('time = {:d} ms'.format( np.int( Time[t] ) ) )
    else:
      time_text.set_text('time = {:d}'.format( np.int( Time[t] ) ) )
    return im,

# Call the animator
print( 'Making movie...' )
anim = animation.FuncAnimation( fig, UpdateFrame, frames = FileArray.shape[0], \
                                interval = 100, blit = True )

anim.save( '{:}_{:}.mp4'.format( ProblemName, Field ), fps = 5 )
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
