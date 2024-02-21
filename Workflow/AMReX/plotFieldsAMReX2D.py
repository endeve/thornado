#!/usr/bin/env python3

import numpy as np
import yt

import GlobalVariables.Settings as gvS

from Utilities.Files import GetFileNumberArray

### Beginning of user input ###

# Specify name of problem

ProblemName = 'YahilCollapse_XCFC'
PlotTitle = ProblemName
suffix = 'png'
figName = 'fig.{:}.{:}'.format( PlotTitle, suffix )
coordinateSystem = 'spherical'

gvS.PlotDirectory = '/Users/nickroberts/thornado/SandBox/AMReX/Applications/YahilCollapse_XCFC/'
#gvS.PlotDirectory = '/home/kkadoogan/Work/Codes/thornado/SandBox/AMReX/Applications/YahilCollapse_XCFC/'

# Specify plot file base name
PlotBaseName = ProblemName + '.plt'

FileNumberArray \
  = GetFileNumberArray \
      ( gvS.PlotDirectory, \
        PlotBaseName, \
        -1, -1, \
        1 )

ds = yt.load( gvS.PlotDirectory \
                + PlotBaseName + str( FileNumberArray[-1] ).zfill( 8 ) )

# Units for axes (dimensionless -> code_length)
LengthUnitX = [ 'code_length', 'code_length' ]

# Set xlabel and ylabel
xLabel = r'$r$'
yLabel = r'$z$'

# Zoom in?
Zoom = 1.0

# Specify field to plot
Field = 'PF_D'

# Label for colorbar (z-axis)
zLabel = r'$\rho\,\left[\mathrm{g\,cm}^{3}\right]$'

# Colormap for plot
cmap = 'viridis'

# Use custom limits for colorbar
UseCustomZmin = False ; zmin = 0.01
UseCustomZmax = False ; zmax = 5.0e1

# Use log scale for z-axis
UseLogScaleZ = True

# Include minor tickmarks (x,Y)
ShowMinorTicksXY = True

# Include minor tickmarks (Z)
ShowMinorTicksZ = True

# Show mesh
ShowMesh = False
MeshAlpha = 0.5

# Overplot contours
OverplotContours = False
nLevels = 5

### End of user input ###

bl = 'boxlib'
blField = (bl,Field)

if ( coordinateSystem == 'spherical' ) :
    axis = 'phi'
else:
    axis = 'z'

slc \
  = yt.SlicePlot \
      ( ds, axis, blField, \
        origin = 'native' )

slc.set_axes_unit( LengthUnitX )

if ( ShowMesh and coordinateSystem != 'spherical' ) :
    slc.annotate_cell_edges \
      ( line_width = 1.0e-12, alpha = MeshAlpha, color = 'black' )

slc.set_cmap( blField, cmap )

slc.set_colorbar_label( Field, zLabel )

if ( not UseCustomZmin ) : zmin = 'min'
if ( not UseCustomZmax ) : zmax = 'max'

slc.set_zlim( Field, zmin = zmin, zmax = zmax )

if ( OverplotContours ) :
    slc.annotate_contour \
      ( blField, levels = nLevels, clim = (zmin,zmax) )

slc.set_log( Field, log = UseLogScaleZ )

slc.set_minorticks         ( Field, ShowMinorTicksXY )
slc.set_colorbar_minorticks( Field, ShowMinorTicksZ  )

slc.set_xlabel( xLabel )
slc.set_ylabel( yLabel )

slc.zoom( Zoom )

slc.save( figName, suffix = suffix, mpl_kwargs = {'bbox_inches':'tight'} )
