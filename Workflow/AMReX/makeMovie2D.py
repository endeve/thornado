#!/usr/bin/env python3

import numpy as np
import yt
from os.path import isfile

from Utilities.Files import GetFileNumberArray

yt.set_log_level( 40 )

### Beginning of user input ###

# Name of problem
ProblemName = 'SedovTaylorBlastWave_XCFC'

PlotTitle = ProblemName

# Directory containing amrex plotfiles

THORNADO_DIR = '/mnt/shared/work/codes/thornado/'

PlotDirectory \
  = THORNADO_DIR + \
      'SandBox/AMReX/Applications/{0:}/{0:}_2D_Cylindrical/nLevels03/' \
      .format( ProblemName )

# Plotfile base name
PlotBaseName = ProblemName + '.plt'

showMesh = True
MeshAlpha = 0.5

xLabel = r'$r$'
yLabel = r'$z$'

# Zoom in?
Zoom = 1.0

# Specify field to plot
Field = 'AF_P'
# Label for colorbar (z-axis)
zLabel = r'$p$'

# Colormap for plot
cmap = 'viridis'

# Use custom limits for colorbar
UseCustomZmin = True ; zmin = 1.0e-10
UseCustomZmax = True ; zmax = 1.0e-1

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

### End of user input ###

# Just used to get bounds
ds = yt.load( PlotDirectory + PlotBaseName + '00000000'  )

# Get lower and higher boundaries and convert them to numpy arrays
xL = ds.domain_left_edge.to_ndarray()
xH = ds.domain_right_edge.to_ndarray()

xL = np.array( [ 0.0, -0.4, 0.0 ] )
xH = np.array( [ 0.5, +0.4, 2.0 * np.pi ] )
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

def createFrame( index, N ):

    fileName = '{:}.png'.format( str( index ).zfill( 3 ) )

    if ( isfile( fileName ) ) : return

    ds = yt.load( PlotDirectory + '{:}{:}' \
                  .format( PlotBaseName, \
                           str( FileNumberArray[index] ).zfill( 8 ) ) )

    slc \
      = yt.SlicePlot \
          ( ds, 'theta', blField, \
            origin = 'native', \
            center = center, \
            width = width )

    slc.set_axes_unit( LengthUnitX )

    slc.annotate_text( [0.7,0.95], \
                       't = {:.2f}'.format( ds.current_time.to_ndarray() ), \
                       coord_system = 'axis', \
                       text_args={'color':'white'})

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
           .format( index / np.float64( N ) * 100.0, fileName ) )

#createFrame( 0 )
N = FileNumberArray.shape[0]
for i in range( N ):
    createFrame( i, N )

'''
To make a movie, generate the *.png files
(be sure to remove the *Render*png files)
and exectute
ffmpeg -f image2 -framerate 30 -i %03d.png -vcodec libx264 -crf 22 mov.2D.mp4
(adapted from https://superuser.com/questions/1462071/ffmpeg-convert-pngs-to-video-files)
'''
