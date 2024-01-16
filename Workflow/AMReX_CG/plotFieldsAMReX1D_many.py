#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )


import GlobalVariables.Settings as gvS
import GlobalVariables.Units    as gvU

from Utilities.GetPlotData  import GetPlotData
from Utilities.Files        import GetFileNumberArray

"""

Default usage, plots last Plotfile in PlotDirectory:

  $ python3 plotFieldsAMReX.py

Alernate usage, plot specific file in PlotDirectory:

  $ python3 plotFieldsAMReX.py 10

  will plot the *00000010 Plotfile

"""

#### ========== User Input ==========

# Specify name of problem
ProblemName = 'YahilCollapse_XCFC'

# Specify title of figure
FigTitle = ProblemName

# Specify directory containing amrex Plotfiles
PlotDirectory = '/Users/nickroberts/thornado/SandBox/AMReX/Applications/YahilCollapse_XCFC/Data_9Lvls_512/'



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
MaxLevel = -1

# Write extra info to screen
Verbose = True

# Use custom limts for y-axis (1D) or colorbar (2D)
UseCustomLimits = False
vmin = 0.0
vmax = 2.0

ShowRefinement = True
RefinementLocations = [ 5.0e+4, 2.5E+4, 1.25E+4, 6.25E+3, 3.125E+3, \
                        1.5625E+3, 7.8125E+2, 3.90625E+2, 1.953125E+2 ]

# Save figure (True) or plot figure (False)
SaveFig = False

# Specify colormap (2D only)
cmap = 'viridis'



FileNumberArray = GetFileNumberArray( PlotDirectory,      \
                                      PlotBaseName        )

nFiles = FileNumberArray.shape[0]

FileNumberSubArray = np.array( [ FileNumberArray[0],        \
                                 FileNumberArray[nFiles//2],\
                                 FileNumberArray[-1]    ]   )




#### ====== End of User Input =======

polar = False
if CoordinateSystem == 'spherical':
    polar = True

ID      = '{:s}_{:s}'.format( ProblemName, Field )
FigName = 'fig.{:s}.png'.format( ID )

# Append "/" to PlotDirectory, if not present
if not PlotDirectory[-1] == '/': PlotDirectory += '/'

gvU.SetSpaceTimeUnits(CoordinateSystem, UsePhysicalUnits)





fig, ax = plt.subplots( 1, 1 )




nFiles = FileNumberSubArray.shape[0]

for f in range(nFiles):

    Data, DataUnit, Time,   \
    X1_C, X2_C, X3_C,       \
    dX1, dX2, dX3, xL, xH = GetPlotData( PlotDirectory,     \
                                         PlotBaseName,      \
                                         Field,             \
                                         argv = ['a',str(FileNumberSubArray[f])],   \
                                         Verbose = Verbose  )



    nX1 = X1_C.shape[0]
    nX2 = X2_C.shape[0]
    nX3 = X3_C.shape[0]


    nDims = 1
    if nX2 > 1: nDims += 1
    if nX3 > 1: nDims += 1

    if nDims == 1:

        ax.plot( X1_C, Data, 'k.' )
        
    else:

        print('Multi-D not implemented yet.')
        
     
     
     
     
if UseLogScale_X:
    ax.set_xscale( 'log' )
    xL = [ max( xL[0], 0.0 + 0.25 * (X1_C[1]-X1_C[0]) ), 0 ]
    
if UseLogScale_Y: ax.set_yscale( 'log' )

if UseCustomLimits: ax.set_ylim( vmin, vmax )

if ShowRefinement:
    bottom, top = plt.ylim()
    ax.plot( (RefinementLocations[:], RefinementLocations[:]), \
             (top, bottom),     \
             scaley = False,    \
             color  = 'red',    \
             zorder = 0,        \
             alpha  = 0.4       )



ax.set_xlim( xL[0], xH[0] )
ax.set_xlabel \
  ( r'$x^{{1}}\ \left[\mathrm{{{:}}}\right]$'.format( gvU.X1Units ), \
    fontsize = 15 )
    
ax.set_ylabel( Field + ' ' + '$'+DataUnit+'$' )
ax.grid()








if SaveFig:

    plt.savefig( FigName, dpi = 300 )

else:

    plt.show()

plt.close()

import os
os.system( 'rm -rf __pycache__ ' )

