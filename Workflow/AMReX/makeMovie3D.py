#!/usr/bin/env python3

import numpy as np
import yt
from os.path import isfile
from yt.visualization.volume_rendering.api import \
    BoxSource, CoordinateVectorSource, GridSource, Scene, create_volume_source

from Utilities.Files import GetFileNumberArray

yt.set_log_level( 40 )

### Beginning of user input ###

# Specify name of problem
ProblemName = 'SedovTaylorBlastWave_XCFC'

Field = 'PF_E'

prefix = '{0:}_3d_{1:}'.format( ProblemName, Field )

# Specify directory containing amrex Plotfiles
THORNADO_DIR = '/mnt/shared/work/codes/thornado/'

PlotDirectory = THORNADO_DIR + 'SandBox/AMReX/Applications/{:}/'.format( ProblemName )

# Specify plot file base name
PlotBaseName = ProblemName + '.plt'

FileNumberArray \
  = GetFileNumberArray \
      ( PlotDirectory, \
        PlotBaseName, \
        -1, -1, \
        1 )

useLogScale = True
def mylog10( x ):
    if ( useLogScale ) :
        return np.log10( x )
    else:
        return x

bounds = ( 1.0e-5, 1.0e-3 )
def linramp( vals, minval, maxval ):
    return ( vals - vals.min() ) / ( vals.max() - vals.min() )

def createFrame( index, xx, overwrite ):

    fileName = 'fig.{0:}_{1:}.png'.format( prefix, str( index ).zfill( 3 ) )

    if ( ( isfile( fileName ) ) & ( not overwrite ) ) : return

    ds = yt.load( PlotDirectory + '{:}{:}' \
                  .format( PlotBaseName, \
                           str( xx[index] ).zfill( 8 ) ) )

    im, sc = yt.volume_render( ds, ('boxlib',Field) )

    cam = sc.get_camera()
    dx, dy, dz = 1.1, 1.1, 1.1
#    dx, dy, dz = 1.0, 1.0, 1.0
    cam.set_width( [ dx, dy, dz ] )
    cam.set_position( [ +1.0, 0.0, 0.0 ], north_vector = [ 0, 0, 1 ] )

    source = sc.get_source()

    # Draw domain boundary
    box_source \
      = BoxSource \
          ( ds.domain_left_edge, \
            ds.domain_right_edge, \
            [ 1.0, 1.0, 1.0, 1.0 ] )
    sc.add_source( box_source )

    # Use custom transfer function
    tf = yt.ColorTransferFunction( mylog10( bounds ) )
    tf.map_to_colormap(
      mylog10( bounds[0] ), \
      mylog10( bounds[1] ), \
      scale = 1.0, \
      colormap = 'YlOrBr', \
      scale_func = linramp )
    source.tfh.set_log( useLogScale )
    source.set_transfer_function( tf )

    # Draw coordinate axes
    coord_source = CoordinateVectorSource()
    sc.add_source( coord_source )

#    # Draw AMR grids
#    grid_source = GridSource( ds.all_data(), alpha = 1.0 )
#    sc.add_source( grid_source )

    sc.save_annotated \
      ( fileName, \
        text_annotate \
          = [ [ (0.05,0.05), \
                't = {:.2f}'.format( ds.current_time ), \
                dict( horizontalalignment = 'left' ) ], \
              [ (0.5,0.95), \
                'Sedov--Taylor Blast Wave', \
                dict( color = 'w', fontsize = '24', \
                      horizontalalignment = 'center') ] ], \
        sigma_clip = 0.7 )

    print( '  [{:.2f}%] Saved {:}' \
           .format( index / np.float64( xx.shape[0] ) * 100.0, fileName ) )

print( '\n  Data directory: {:}'.format( PlotDirectory ) )
xx = FileNumberArray#[0:-1:10]
overwrite = True
for i in range( 0, FileNumberArray.shape[0], 1 ):
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
ffmpeg -loglevel quiet -f image2 -framerate 30 -i %03d.png -vcodec libx264 -crf 22 mov.3D.mp4
(adapted from https://superuser.com/questions/1462071/ffmpeg-convert-pngs-to-video-files)
'''
