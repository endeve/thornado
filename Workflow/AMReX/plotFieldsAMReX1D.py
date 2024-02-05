#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

import GlobalVariables.Settings as gvS
import GlobalVariables.Units    as gvU

from Utilities.GetPlotData              import GetPlotData
from Utilities.RefinementBoundaryFinder import FindRefinementBoundaries

"""

Default usage, plots last Plotfile in PlotDirectory:

  $ python3 plotFieldsAMReX.py

Alernate usage, plot specific file in PlotDirectory:

  $ python3 plotFieldsAMReX.py 10

  will plot the *00000010 Plotfile

"""

#### ========== User Input ==========

DataType = 'AMReX'

# Specify name of problem
ProblemName = 'AdiabaticCollapse_XCFC'

# Specify title of figure
FigTitle = '{:}'.format( ProblemName )

# Specify directory containing amrex Plotfiles
PlotDirectory = '/Users/nickroberts/thornado/SandBox/AMReX/Applications/AdiabaticCollapse_XCFC/'
#PlotDirectory \
#  = '/home/kkadoogan/Work/Codes/thornado/\
#SandBox/AMReX/dgExperiments_Euler_Relativistic_IDEAL/'

# Specify plot file base name
PlotBaseName = ProblemName + '.plt'

# Specify field to plot
Field = 'PF_D'

# Specify to plot in log-scale
UseLogScale_X  = True
UseLogScale_Y  = True

# Specify whether or not to use physical units
UsePhysicalUnits = True

# Specify coordinate system (currently supports 'cartesian' and 'spherical')
CoordinateSystem = 'spherical'

# Max level of refinement to plot (-1 plots leaf elements)
MaxLevel = -1

# Write extra info to screen
Verbose = True

# Use custom limts for y-axis
UseCustomLimits = False
vmin = 0.0
vmax = 2.0

ShowRefinement = True

# Save figure (True) or plot figure (False)
SaveFig = False

#### ====== End of User Input =======

ID      = '{:s}_{:s}'.format( ProblemName, Field )
FigName = 'fig.{:s}.png'.format( ID )

# Append "/" to PlotDirectory, if not present
if not PlotDirectory[-1] == '/': PlotDirectory += '/'

gvU.SetSpaceTimeUnits( CoordinateSystem, UsePhysicalUnits )

Data, DataUnit, Time, X1_C, X2_C, X3_C, dX1, dX2, dX3 \
  = GetPlotData \
      ( PlotDirectory      , \
        PlotBaseName       , \
        Field              , \
        argv = argv        , \
        DataType = DataType, \
        Verbose = Verbose )

xL = X1_C[0 ] - 0.5 * dX1[0 ]
xH = X1_C[-1] + 0.5 * dX1[-1]

nX1 = X1_C.shape[0]
nX2 = X2_C.shape[0]
nX3 = X3_C.shape[0]

RefinementLocations = FindRefinementBoundaries( dX1 )

RefinementLocations = []
x = xL
for i in range( nX1 ):
    x += dX1[i]
    RefinementLocations.append( x )
xRef = np.array( RefinementLocations )

nDims = 1
if nX2 > 1: nDims += 1
if nX3 > 1: nDims += 1

fig, ax = plt.subplots( 1, 1 )

ax.plot( X1_C, Data, 'k.' )

if UseLogScale_X:
    ax.set_xscale( 'log' )
    xL = [ max( xL, 0.0 + 0.25 * dX1[0] ), 0 ]

if UseLogScale_Y: ax.set_yscale( 'log' )

if UseCustomLimits: ax.set_ylim( vmin, vmax )

for xx in xRef:
    ax.axvline( xx )
try:
    ax.set_xlim( xLim )
except:
    ax.set_xlim( xL, xH )

try:
    ax.set_xlabel( xLabel, fontsize = 15 )
except:
    ax.set_xlabel \
      ( r'$x\ \left[\mathrm{{{:}}}\right]$'.format( gvU.X1Units ), \
        fontsize = 15 )

try:
    ax.set_ylabel( yLabel, fontsize = 15 )
except:
    ax.set_ylabel( Field + ' ' + '$'+DataUnit+'$' )

ax.grid()

try:
    ax.set_title( figTitle )
except:
    ax.set_title( r'$\texttt{{{:}}}$'.format( FigTitle ), fontsize = 15 )

if SaveFig:

    plt.savefig( FigName, dpi = 300 )

else:

    plt.show()

plt.close()

import os
os.system( 'rm -rf __pycache__ ' )
