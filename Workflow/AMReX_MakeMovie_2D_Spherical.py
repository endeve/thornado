#!/usr/local/bin/python3

"""
Given a directory with 2D AMReX plot-files
in spherical-polar coordinates, create a data
file containing all the values of a specified
variable for each time-step and use that to
generate a movie.
"""

import yt
import numpy as np
import subprocess
from os import listdir
from os.path import isfile
from sys import argv, exit
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib.ticker as ticker

yt.funcs.mylog.setLevel(0) # Suppress initial yt output to screen

# --- Get user's HOME directory ---
HOME = subprocess.check_output( ["echo $HOME"], shell = True)
HOME = HOME[:-1].decode( "utf-8" ) + '/'

############# User input #############
DataDirectory = HOME + 'Desktop/'

ProblemName = 'SASI' # Only used for name of movie file

# Create new variables following method used near line 98
VariableToPlot = 'Entropy'

DataFileName     = 'MovieData_{:}.dat'.format( VariableToPlot )
TimeFileName     = 'MovieTime.dat'
UseLogScale      = False      # Do you want your movie in log scale?
UsePhysicalUnits = True      # Are you using physical units?
Relativistic     = True      # Are you plotting results from relativistic hydro?
cmap             = 'Purples' # Color scheme for movie
UseCustomTicks   = True      # Define limits for colorbar near line 166
zAxisVertical    = True     # Orient z-axis
############# End of user input #############

if( UsePhysicalUnits ):
    c = 2.99792458e10
    Centimeter = 1.0e5 # Centimeters per kilometer
else:
    c = 1.0
    Centimeter = 1.0

# Get last plotfile in directory
FileArray = np.sort(np.array( [ file for file in listdir( DataDirectory ) ] ) )
FileList = []
for iFile in range( FileArray.shape[0] ):
    sFile = FileArray[iFile]
    if( sFile[0:8] == 'thornado' ):
        FileList.append( sFile )
FileArray = np.array( FileList )
File = FileArray[-1]

ds = yt.load( '{:}'.format( DataDirectory + File ) )

MaxLevel = ds.index.max_level
Time     = ds.current_time
nX       = ds.domain_dimensions
xL       = ds.domain_left_edge
xH       = ds.domain_right_edge

Overwrite = True
if( isfile( DataFileName ) ):
    Overwrite = input( 'File: "{:}" exists. overwrite? (Y/N): '.format \
                  ( DataFileName ) )
    if( not Overwrite == 'Y' ):
        print( 'Not overwriting file, using file {:} for movie.'.format( \
          DataFileName ) )
        Overwrite = False
    else:
        Overwrite = True

if( Overwrite ):
    # Put all time-slices into one array to use for movie making
    Data = np.empty( (FileArray.shape[0],nX[0],nX[1]), float )
    Time = np.empty( FileArray.shape[0], float )
    print( 'Generating data file: {:}...'.format( DataFileName ) )
    for i in range( FileArray.shape[0] ):
        print( '{:}/{:}'.format( i+1, FileArray.shape[0] ) )
        ds = yt.load( '{:}'.format( DataDirectory + FileArray[i] ) )

        CoveringGrid \
          = ds.covering_grid \
              ( level           = MaxLevel, \
                left_edge       = xL, \
                dims            = nX * 2**MaxLevel, \
                num_ghost_zones = nX[0] )

        if( VariableToPlot == 'Entropy' ):
            PF_D  = CoveringGrid['PF_D' ].to_ndarray()[:,:,0]
            AF_P  = CoveringGrid['AF_P' ].to_ndarray()[:,:,0]
            AF_Gm = CoveringGrid['AF_Gm'].to_ndarray()[:,:,0]
            Data[i] = AF_P / PF_D**AF_Gm
        elif( VariableToPlot == 'LorentzFactor' ):
            PF_V1 = CoveringGrid['PF_V1'   ].to_ndarray()[:,:,0] * Centimeter
            PF_V2 = CoveringGrid['PF_V2'   ].to_ndarray()[:,:,0] * Centimeter
            PF_V3 = CoveringGrid['PF_V3'   ].to_ndarray()[:,:,0] * Centimeter
            GF_g1 = CoveringGrid['GF_Gm_11'].to_ndarray()[:,:,0]
            GF_g2 = CoveringGrid['GF_Gm_22'].to_ndarray()[:,:,0]
            GF_g3 = CoveringGrid['GF_Gm_33'].to_ndarray()[:,:,0]
            Data[i] \
              = 1 / np.sqrt( 1.0 - ( GF_g1 * PF_V1**2 \
                                     + GF_g2 * PF_V2**2 \
                                     + GF_g3 * PF_V3**2 ) / c**2 )
        elif( VariableToPlot == 'SpecificEnthalpy' ):
            PF_D = CoveringGrid['PF_D'].to_ndarray()[:,:,0]
            PF_E = CoveringGrid['PF_E'].to_ndarray()[:,:,0]
            AF_P = CoveringGrid['AF_P'].to_ndarray()[:,:,0]
            if( Relativistic ): # Plot h/c^2
                Data[i] = ( c**2 + ( PF_E + AF_P ) / PF_D ) / c**2
            else:
                Data[i] = ( PF_E + AF_P ) / PF_D
        else:
            Data[i] = CoveringGrid[VariableToPlot].to_ndarray()[:,:,0]
            if( VariableToPlot[-2] == 'V' ):
                Data[i] *= Centimeter

        Time[i] = ds.current_time

    # Save multi-D array with np.savetxt. Taken from:
    # https://stackoverflow.com/questions/3685265/how-to-write-a-multidimensional-array-to-a-text-file

    with open( DataFileName, 'w' ) as FileOut:
        FileOut.write( '# Array shape: {0}\n'.format( Data.shape ) )

        # Iterating through an n-dimensional array produces slices along
        # the last axis. This is equivalent to Data[i] in this case
        for TimeSlice in Data:
            np.savetxt( FileOut, TimeSlice )
            FileOut.write( '# New slice\n' )

    np.savetxt( TimeFileName, Time )

print( 'Reading in data file: {:}...'.format( DataFileName ) )
Data = np.loadtxt( DataFileName ).reshape( \
                                    (FileArray.shape[0],nX[0],nX[1]) )
Time = np.loadtxt( TimeFileName )

fig = plt.figure()
ax  = fig.add_subplot( 111, polar = True )
xL  = xL.to_ndarray()
xH  = xH.to_ndarray()
X1  = np.linspace( xL[0], xH[0], nX[0] )
X2  = np.linspace( xL[1], xH[1], nX[1] )
theta, r = np.meshgrid( X2, X1 )

if( UseLogScale ):
    from matplotlib.colors import LogNorm
    norm = LogNorm()
    def f(t):
        return np.abs( Data[t] )
else:
    norm = None
    def f(t):
        return Data[t]

if( UseCustomTicks ):
    vmin = 1.0
    vmax = max( -np.inf, np.max( Data ) )
    if( UseLogScale ):
        ticks = np.logspace( np.log10( vmin ), np.log10( vmax ), 5 )
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
im = ax.pcolormesh( theta, r, f(0)[:-1,:-1], \
                    cmap = cmap, \
                    vmin = vmin, vmax = vmax, \
                    norm = norm )
ax.set_thetamin( 180.0/np.pi * X2[0 ] )
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
    fig.colorbar( im )

def UpdateFrame(t):
    im.set_array( f(t)[:-1,:-1].flatten() )
    if( UsePhysicalUnits ):
      time_text.set_text('time = {:d} ms'.format( np.int( Time[t] ) ) )
    else:
      time_text.set_text('time = {:d}'.format( np.int( Time[t] ) ) )
    return im,

# Call the animator
print( 'Making movie...' )
anim = animation.FuncAnimation( fig, UpdateFrame, frames = FileArray.shape[0], \
                                interval = 100, blit = True )

anim.save( '{:}_{:}.mp4'.format( ProblemName, VariableToPlot ), fps = 5 )
plt.close()
