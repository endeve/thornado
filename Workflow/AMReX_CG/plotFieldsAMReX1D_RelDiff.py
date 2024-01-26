#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )


import GlobalVariables.Settings as gvS
import GlobalVariables.Units    as gvU

from Utilities.GetPlotData  import GetPlotData
from Utilities.GetFrameData import GetFrameData

"""

Default usage, plots last Plotfile in PlotDirectory:

  $ python3 plotFieldsAMReX.py

Alernate usage, plot specific file in PlotDirectory:

  $ python3 plotFieldsAMReX.py 10

  will plot the *00000010 Plotfile

"""

#### ========== User Input ==========

# Specify name of problem
ProblemName = 'AdiabaticCollapse_XCFC'

# Specify title of figure
FigTitle = ProblemName

# Specify directory containing amrex Plotfiles
gvS.nDirs = 2

PlotDirectories = ['None']*gvS.nDirs

#PlotDirectories[0] = '/Users/nickroberts/thornado/SandBox/AMReX/Applications/AdiabaticCollapse_XCFC'
#PlotDirectories[1] = '/Users/nickroberts/thornado/SandBox/AdiabaticCollapse_XCFC/Output'

PlotDirectories[0] = '/Users/nickroberts/thornado_clean/thornado/SandBox/AdiabaticCollapse_XCFC/Output'
PlotDirectories[1] = '/Users/nickroberts/thornado_clean/thornado/SandBox/AMReX/Applications/AdiabaticCollapse_XCFC'
    
gvS.DataType = ['None']*gvS.nDirs
#gvS.DataType[0] = 'AMReX'
#gvS.DataType[1] = 'Native'
 
gvS.DataType[0] = 'Native'
gvS.DataType[1] = 'AMReX'
 
FrameNumber = 0 # -1 -> argv, or last frame



# Specify plot file base name
PlotBaseName = ProblemName + '.plt'

# Specify field to plot
Field = 'PF_D'

# Specify to plot in log-scale
UseLogScale_X  = True
UseLogScale_Y  = True
UseLogScale_2D = False

# Specify whether or not to use physical units
UsePhysicalUnits = False

# Specify coordinate system (currently supports 'cartesian' and 'spherical')
CoordinateSystem = 'cartesian'

# Max level of refinement to plot (-1 plots leaf elements)
gvS.MaxLevel = -1

# Include initial conditions in movie?
gvS.ShowIC = True

gvS.PlotMesh = False

# Write extra info to screen
gvS.Verbose = True

# Use custom limts for y-axis (1D) or colorbar (2D)
gvS.UseCustomLimits = True
gvS.vmin = 1.0e-16
gvS.vmax = 1.0e1

gvS.MovieRunTime = 10.0 # seconds

gvS.ShowRefinement = True
gvS.RefinementLevels = 7

gvS.ReferenceBounce = False


gvS.amr = True
    
    
# Save figure (True) or plot figure (False)
SaveFig = False

# Specify colormap (2D only)
cmap = 'viridis'







#### ====== End of User Input =======

polar = False
if CoordinateSystem == 'spherical':
    polar = True

DataDirectories = ['None']*gvS.nDirs

DataDirectories[0] = 'DataDirectories/{:s}_DWN'.format( ProblemName )
DataDirectories[1] = 'DataDirectories/{:s}_DWNB'.format( ProblemName )

ID            = '{:s}_{:s}'.format( ProblemName, Field )
gvS.MovieName     = 'mov.{:s}_RelDiff.mp4'.format( ID )

# Append "/" to PlotDirectory, if not present
for i in range(gvS.nDirs):
    if not PlotDirectories[i][-1] == '/': PlotDirectories[i] += '/'


gvU.SetSpaceTimeUnits(CoordinateSystem, UsePhysicalUnits)


Data = ['None']*gvS.nDirs
DataUnits = ['None']*gvS.nDirs
X1 = ['None']*gvS.nDirs
X2 = ['None']*gvS.nDirs
X3 = ['None']*gvS.nDirs
Time = ['None']*gvS.nDirs
dX1 = ['None']*gvS.nDirs
dX2 = ['None']*gvS.nDirs
dX3 = ['None']*gvS.nDirs
xL = ['None']*gvS.nDirs
xH = ['None']*gvS.nDirs


for i in range(gvS.nDirs):

    Data[i], DataUnits[i], Time[i],    \
    X1[i], X2[i], X3[i],              \
    dX1[i], dX2[i], dX3[i] = GetPlotData( PlotDirectories[i],     \
                                                        PlotBaseName,           \
                                                        Field,                  \
                                                        gvS.DataType[i],        \
                                                        FrameNumber,            \
                                                        argv,                   \
                                                        gvS.Verbose             )


    xL[i] = X1[i][0 ] - 0.5 * dX1[i][0 ]
    xH[i] = X1[i][-1] + 0.5 * dX1[i][-1]



if (any(X1[0] != X1[1])):
    msg = "\nRaidal meshes do not match. \n"
    msg +=" Can not do a proper relative difference.\n"
    exit(msg)


nX1 = X1[0].shape[0]
nX2 = X2[0].shape[0]
nX3 = X3[0].shape[0]


nDims = 1
if nX2 > 1: nDims += 1
if nX3 > 1: nDims += 1

if nDims == 1:

    fig, ax = plt.subplots( 1, 1 )
    
    
    RelDiff = abs(Data[0]-Data[1])/abs(Data[0])

    ax.plot( X1[0], RelDiff, 'k.' )
    
    
    if UseLogScale_X:
        ax.set_xscale( 'log' )
        print(xL[0][0], 0.0 + 0.25 * dX1[0][0])
        xLA = [ 0.0 + 0.25 * dX1[0][0], 0 ]
        
    if UseLogScale_Y: ax.set_yscale( 'log' )
    
    if UseCustomLimits: ax.set_ylim( vmin, vmax )
    
#    if ShowRefinement:
#        bottom, top = plt.ylim()
#        ax.plot( (RefinementLocations[:], RefinementLocations[:]), \
#                 (top, bottom),     \
#                 scaley = False,    \
#                 color  = 'red',    \
#                 zorder = 0,        \
#                 alpha  = 0.4       )



    ax.set_xlim( xL[0][0], xH[0][0] )
    ax.set_xlabel \
      ( r'$x^{{1}}\ \left[\mathrm{{{:}}}\right]$'.format( gvU.X1Units ), \
        fontsize = 15 )
        
    ax.set_ylabel( Field + ' ' + '$'+DataUnitA+'$' )
    ax.grid()

else:

    print('Multi-D not implemented yet.')






if SaveFig:

    plt.savefig( FigName, dpi = 300 )

else:

    plt.show()

plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
