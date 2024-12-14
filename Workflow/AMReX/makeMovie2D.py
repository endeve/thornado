#!/usr/bin/env python3

import numpy as np
import yt
from os.path import isfile
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

from Utilities.Files import GetFileNumberArray

yt.set_log_level( 40 )

### Beginning of user input ###

# Name of problem
ProblemName = 'SedovTaylorBlastWave_XCFC'

PlotTitle = ProblemName

# Directory containing amrex plotfiles

THORNADO_DIR = '/mnt/shared/work/codes/thornado/'

suffix = ''
PlotDirectory \
  = THORNADO_DIR + \
      'SandBox/AMReX/Applications/{0:}/' \
      .format( ProblemName, suffix )

# Plotfile base name
PlotBaseName = ProblemName + '.plt'

showMesh = True
MeshAlpha = 0.5

xLabel = r'$r$'
yLabel = r'$z$'

# Zoom in?
Zoom = 1.0

# Specify field to plot
Field = 'PF_D'
# Label for colorbar (z-axis)
zLabel = r'$\rho$'

# Colormap for plot
cmap = 'viridis'

# Use custom limits for colorbar
UseCustomZmin = True ; zmin = +1.0e+0
UseCustomZmax = True ; zmax = +1.0e+3

# Use log scale for z-axis
UseLogScaleZ = True

# Include minor tickmarks (x,Y)
ShowMinorTicksXY = True

# Include minor tickmarks (Z)
ShowMinorTicksZ = True

# Show mesh
ShowMesh  = True
MeshAlpha = 0.5

# Overplot contours
OverplotContours = False
nLevels = 5

prefix = '{0:}_2d_{1:}{2:}'.format( ProblemName, Field, suffix )

### End of user input ###

# Just used to get bounds
ds = yt.load( PlotDirectory + PlotBaseName + '00000000'  )

# Get lower and higher boundaries and convert them to numpy arrays
xL = ds.domain_left_edge.to_ndarray()
xH = ds.domain_right_edge.to_ndarray()

#xL = np.array( [ 0.0, -0.4, 0.0 ] )
#xH = np.array( [ 0.5, +0.4, 2.0 * np.pi ] )
# Set center and width of plot window
center = 0.5 * ( xL + xH )
width  = 1.0 * ( xH - xL )

LengthUnitX = [ 'code_length', 'code_length' ]

FileNumberArray \
  = GetFileNumberArray \
      ( PlotDirectory, \
        PlotBaseName, \
        -1, -1, 1 )

blField = ('boxlib',Field)

from yt.visualization.api import PlotCallback
class timeLabel( PlotCallback ):
    _type_name = "timeLabel"

    def __init__( self, text, x, y ):
        self.text = text
        self.position = ( x, y )

    def __call__( self, plot ):
        plot._axes.annotate(
            self.text,
            xy = self.position,
            xycoords = "figure fraction",
            bbox = { 'facecolor' : 'white', 'pad' : 2, 'alpha' : 0.6 }, \
            fontsize = 16,
            color = 'black' )

def createFrame( index, xx, overwrite ):

    fileName = 'fig.{0:}_{1:}.png'.format( prefix, str( index ).zfill( 3 ) )

    if ( ( isfile( fileName ) ) & ( not overwrite ) ) : return

    ds = yt.load( PlotDirectory + '{:}{:}' \
                  .format( PlotBaseName, \
                           str( xx[index] ).zfill( 8 ) ) )

    slc \
      = yt.SlicePlot \
          ( ds, 'theta', blField, \
            origin = 'native', \
            center = center, \
            width = width )

    slc.set_axes_unit( LengthUnitX )
    slc.annotate_timeLabel \
      ( r'$t = {:.2f}$'.format( np.float64( ds.current_time.to_ndarray() ) ), \
        0.55, 0.9 )

    #slc.annotate_grids()

    if ( showMesh ) :
        slc.annotate_cell_edges \
            ( line_width = 1.0e-12, alpha = MeshAlpha, color = 'black' )

    slc.set_cmap( blField, cmap )

    slc.set_colorbar_label( Field, zLabel )

    global zmin, zmax

    if ( not UseCustomZmin ) : zmin = 'min'
    if ( not UseCustomZmax ) : zmax = 'max'

    slc.set_zlim( Field, zmin = zmin, zmax = zmax )

    if ( OverplotContours ) :
        slc.annotate_contour \
          ( blField, levels = nLevels, clim = (zmin,zmax) )

    slc.set_log( Field, log = UseLogScaleZ )

    slc.set_minorticks         ( Field, ShowMinorTicksXY )
    slc.set_colorbar_minorticks( Field, ShowMinorTicksZ  )

    slc.annotate_title( PlotTitle )

    slc.set_xlabel( xLabel )
    slc.set_ylabel( yLabel )

    slc.zoom( Zoom )

    slc.save( fileName )

    print( '  [{:.2f}%] Saved {:}' \
           .format( index / np.float64( xx.shape[0] ) * 100.0, fileName ) )

print( '\n  Data directory: {:}'.format( PlotDirectory ) )
xx = FileNumberArray#[0:-1:10]
overwrite = False
for i in range( xx.shape[0] ):
    createFrame( i, xx, overwrite )

movName = 'mov.{:}.mp4'.format( prefix )
import os
os.system( \
'ffmpeg -loglevel quiet -f image2 -framerate 30 -i fig.{0:}_%03d.png -vcodec libx264 -crf 22 {1:}'.format( prefix, movName ) )
print( '\n  Saved {:}'.format( movName ) )

'''
To make a movie, generate the *.png files
(be sure to remove the *Render*png files)
and exectute
ffmpeg -loglevel quiet -f image2 -framerate 30 -i %03d.png -vcodec libx264 -crf 22 mov.2D.mp4
(adapted from https://superuser.com/questions/1462071/ffmpeg-convert-pngs-to-video-files)
'''
