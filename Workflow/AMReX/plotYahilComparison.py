#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

import GlobalVariables.Settings as gvS
import GlobalVariables.Units    as gvU

from Utilities.YahilProfile import LoadYahilProfile, GetYahilValues
from Utilities.GetPlotData import GetPlotData


time = 51.0 # s

gamma = 1.3
central_density = 7e9 # g cm^-3
pressure = 6e27 # erg cm^-3

kappa = pressure/central_density**gamma


YahilFile = '/Users/nickroberts/thornado/Workflow/AMReX_MR/YahilHomologousCollapse_Gm_130.dat'

# Specify name of problem
ProblemName = 'YahilCollapse_XCFC'

# Specify title of figure
FigTitle = ProblemName

# Specify directory containing amrex Plotfiles
PlotDirectory = '/Users/nickroberts/thornado/SandBox/AMReX/Applications/YahilCollapse_XCFC/Data_9Lvls_512/'

# Specify plot file base name
PlotBaseName = ProblemName + '.plt'

# Specify to plot in log-scale
UseLogScale_X  = True
UseLogScale_Y  = True
UseLogScale_2D = False

# Specify whether or not to use physical units
UsePhysicalUnits = False

# Specify coordinate system (currently supports 'cartesian' and 'spherical')
CoordinateSystem = 'cartesian'

# Max level of refinement to plot (-1 plots leaf elements)
MaxLevel = -1

# Write extra info to screen
Verbose = True

# Use custom limts for y-axis (1D) or colorbar (2D)
UseCustomLimits = False
vmin = 0.0
vmax = 2.0

# Save figure (True) or plot figure (False)
SaveFig = False

# Specify colormap (2D only)
cmap = 'viridis'






#### ====== End of User Input =======

polar = False
if CoordinateSystem == 'spherical':
    polar = True

# Specify field to plot
FieldA = 'PF_D'
FieldB = 'PF_V1'

FigName = 'fig.YahilTest.png'

# Append "/" to PlotDirectory, if not present
if not PlotDirectory[-1] == '/': PlotDirectory += '/'

gvU.SetSpaceTimeUnits(CoordinateSystem, UsePhysicalUnits)


ThorDen, ThorDenUnits, Time, \
X1_C, X2_C, X3_C,               \
dX1, dX2, dX3, xL, xH = GetPlotData( PlotDirectory,     \
                                     PlotBaseName,      \
                                     'PF_D',            \
                                     argv = argv,       \
                                     Verbose = Verbose  )

ThorVel, ThorVelUnits, Time, \
X1_C, X2_C, X3_C,               \
dX1, dX2, dX3, xL, xH = GetPlotData( PlotDirectory,     \
                                     PlotBaseName,      \
                                     'PF_V1',            \
                                     argv = argv,       \
                                     Verbose = Verbose  )



YahilProfile = LoadYahilProfile( YahilFile )

Yahil_Density, Yahil_Velocity = GetYahilValues( X1_C, kappa, gamma, Time, YahilProfile)


nX = ['None']*3
# Re-define nX
nX[0] = X1_C.shape[0]
nX[1] = X2_C.shape[0]
nX[2] = X2_C.shape[0]

nDims = 1
if nX[1] > 1: nDims += 1
if nX[2] > 1: nDims += 1


if nDims == 1:
    
    
    fig, ax = plt.subplots( 2, 1 )

    ax[0].plot( X1_C, ThorDen, 'r-',label=r'Thornado' )
    ax[0].plot( X1_C, Yahil_Density, 'b--',label=r'Profile')
    ax[1].plot( X1_C, ThorVel, 'r-',label=r'Thornado' )
    ax[1].plot( X1_C, Yahil_Velocity, 'b--',label=r'Profile')
    
    if UseLogScale_X:
        ax[0].set_xscale( 'log' )
        ax[1].set_xscale( 'log' )
        xL = [ max( xL[0], 0.0 + 0.25 * dX1[0] ), 0 ]
    if UseLogScale_Y:
        ax[0].set_yscale( 'log' )
#        ax[1].set_yscale( 'log' )
    if UseCustomLimits:
        ax[0].set_ylim( vmin, vmax )
        ax[1].set_ylim( vmin, vmax )
    ax[0].set_xlim( xL[0], xH[0] )
    ax[1].set_xlim( xL[0], xH[0] )
    ax[0].set_xlabel \
      ( r'$x^{{1}}\ \left[\mathrm{{{:}}}\right]$'.format( gvU.X1Units ), \
        fontsize = 15 )
    ax[1].set_xlabel \
      ( r'$x^{{1}}\ \left[\mathrm{{{:}}}\right]$'.format( gvU.X1Units ), \
        fontsize = 15 )
#    ax[0].set_ylabel( FieldA + ' ' + ThorDenUnits )
#    ax[1].set_ylabel( FieldB + ' ' + ThorVelUnits )
    
    ax[0].legend()

    ax[0].grid()
    ax[1].grid()


ax[0].set_title( r'$\texttt{{{:}}}$'.format( FigTitle ) )

if SaveFig:

    plt.savefig( FigName, dpi = 300 )

else:

    plt.show()

plt.close()

import os
os.system( 'rm -rf __pycache__ ' )




