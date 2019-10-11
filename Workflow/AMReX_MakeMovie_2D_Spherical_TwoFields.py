#!/usr/local/bin/python3

"""
Generate a movie showing two fields side-by-side
from data files created in MakeDataFiles.py.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

############# User input #############
ProblemName = 'SASI' # Only used for name of movie file

# Define custom variables as seen near line 40
Field1 = 'Entropy'
Field2 = 'PF_D'

DataFileName1 = 'MovieData_{:}.dat'.format( Field1 )
DataFileName2 = 'MovieData_{:}.dat'.format( Field2 )

UseLogScale1 = False  # Do you want your first field plotted in log scale?
UseLogScale2 = True  # Do you want your second field plotted in log scale?

cmap1 = 'Purples' # Color scheme for first movie
cmap2 = 'jet'     # Color scheme for second movie

UseCustomTicks = True # Define limits for colorbar near line 166

zAxisVertical = True # Orient z-axis
############# End of user input #############

from MakeDataFiles import *

MakeDataFile( Field1, DataFileName1 )
MakeDataFile( Field2, DataFileName2 )

print( 'Reading in data files...' )
Data1 = np.loadtxt( DataFileName1 ).reshape( \
                                    (FileArray.shape[0],nX[0],nX[1]) )
Data2 = np.loadtxt( DataFileName2 ).reshape( \
                                    (FileArray.shape[0],nX[0],nX[1]) )
Time = np.loadtxt( TimeFileName )

fig = plt.figure( figsize = (8,6) )
ax  = fig.add_subplot( 111, polar = True )
xL  = xL.to_ndarray()
xH  = xH.to_ndarray()
X1  = np.linspace( xL[0], xH[0], nX[0] )
X2  = np.linspace( xL[1], xH[1], nX[1] )
theta1, r = np.meshgrid( X2, X1 )
theta2, r = np.meshgrid( 2.0 * np.pi - X2, X1 )

if( UseLogScale1 ):
    from matplotlib.colors import LogNorm
    norm1 = LogNorm()
    def f1(t):
        return np.abs( Data1[t] )
else:
    norm1 = None
    def f1(t):
        return Data1[t]

if( UseLogScale2 ):
    from matplotlib.colors import LogNorm
    norm2 = LogNorm()
    def f2(t):
        return np.abs( Data2[t] )
else:
    norm2 = None
    def f2(t):
        return Data2[t]

if( UseCustomTicks ):

    vmin1 = min( +np.inf, np.min( Data1 ) )
    vmax1 = max( -np.inf, np.max( Data1 ) )
    vmin2 = min( +np.inf, np.min( Data2 ) )
    vmax2 = max( -np.inf, np.max( Data2 ) )
    if( UseLogScale1 ):
        ticks1 = np.logspace( np.log10( vmin1 ), np.log10( vmax1 ), 5 )
    else:
        ticks1 = np.linspace( vmin1, vmax1, 5 )

    if( UseLogScale2 ):
        ticks2 = np.logspace( np.log10( vmin2 ), np.log10( vmax2 ), 5 )
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

# Taken from:
# https://brushingupscience.com/2016/06/21/matplotlib-animations-the-easy-way/
im1 = ax.pcolormesh( theta1, r, f1(0)[:-1,:-1], \
                     cmap = cmap1, \
                     vmin = vmin1, vmax = vmax1, \
                     norm = norm1 )
im2 = ax.pcolormesh( theta2, r, f2(0)[:-1,:-1], \
                     cmap = cmap2, \
                     vmin = vmin2, vmax = vmax2, \
                     norm = norm2 )
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

def UpdateFrame(t):
    im1.set_array( f1(t)[:-1,:-1].flatten() )
    im2.set_array( f2(t)[:-1,:-1].flatten() )
    if( UsePhysicalUnits ):
      time_text.set_text('time = {:d} ms'.format( np.int( Time[t] ) ) )
    else:
      time_text.set_text('time = {:d}'.format( np.int( Time[t] ) ) )
    return im1, im2,

# Call the animator
print( 'Making movie...' )
anim = animation.FuncAnimation( fig, UpdateFrame, frames = FileArray.shape[0], \
                                interval = 100, blit = True )

anim.save( '{:}_movie.mp4'.format( ProblemName ), fps = 5 )
plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
