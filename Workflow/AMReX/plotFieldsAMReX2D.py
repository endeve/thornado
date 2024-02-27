#!/usr/bin/env python3

import numpy as np
import yt

from Utilities.Files import GetFileNumberArray

### Beginning of user input ###

# Specify name of problem
ProblemName = 'RiemannProblem2D_dZB2002'

# Specify directory containing amrex Plotfiles
#PlotDirectory = '/Users/nickroberts/thornado/SandBox/AMReX/Applications/YahilCollapse_XCFC/'
PlotDirectory = '/home/kkadoogan/Work/Codes/thornado/SandBox/AMReX/dgExperiments_Euler_Relativistic_IDEAL/'

# Specify plot file base name
PlotBaseName = ProblemName + '.plt'

FileNumberArray \
  = GetFileNumberArray \
      ( PlotDirectory, \
        PlotBaseName, \
        -1, -1, \
        1 )

# For a single file, replace `FileNumberArray[-1]` with cycle number,
# e.g., plt00001310 --> str( 1310 )
ds = yt.load( PlotDirectory + PlotBaseName + FileNumberArray[-1].zfill( 8 ) )

# Get lower and higher boundaries and convert them to numpy arrays
xL = ds.domain_left_edge.to_ndarray()
xH = ds.domain_right_edge.to_ndarray()

'''
Get list of plotting options here:
https://yt-project.org/doc/reference/api/yt.visualization.plot_container.html
'''

# Specify field to plot
Field = 'AF_P'

# Colormap for plot
cmap = 'Purples'

# Label for colorbar
cLabel = r'$p$'

# Use custom limits for x-axis
UseCustomZmin = True ; zmin = 0.01
UseCustomZmax = True ; zmax = 5.0e1

UseLogScale_Z = True

# Set xlabel and ylabel
xLabel = r'$x$'
yLabel = r'$y$'

### End of user input ###

bl = 'boxlib'
blField = (bl,Field)

slc = yt.SlicePlot( ds, 'z', blField, origin = 'native' )

slc.set_cmap( blField, cmap )

slc.set_colorbar_label( Field, cLabel )

if ( not UseCustomZmin ) : zmin = 'min'
if ( not UseCustomZmax ) : zmax = 'max'

slc.set_zlim( Field, zmin = zmin, zmax = zmax )

#slc.annotate_contour( blField, levels = 5, clim = (zmin,zmax) )

if ( UseLogScale_Z ) : slc.set_log( Field, log = UseLogScale_Z )

slc.set_minorticks         ( Field, True )
slc.set_colorbar_minorticks( Field, True )

slc.set_xlabel( xLabel )
slc.set_ylabel( yLabel )

#slc._setup_plots()
#plt = slc.plots[bl,Field]
#ax = plt.axes
#ax.set_xlim( xL[0], xH[0] )
#ax.set_ylim( xL[1], xH[1] )

slc.save()
