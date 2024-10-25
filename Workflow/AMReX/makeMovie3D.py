#!/usr/bin/env python3

import numpy as np
import yt
from yt.visualization.volume_rendering.api import \
    BoxSource, CoordinateVectorSource, GridSource, Scene, create_volume_source

from Utilities.Files import GetFileNumberArray

yt.set_log_level( 40 )

### Beginning of user input ###

# Specify name of problem
ProblemName = 'SedovTaylorBlastWave_XCFC'

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

bounds = ( 1.0e-5, 1.0e-2 )
def linramp( vals, minval, maxval ):
    return ( vals - vals.min() ) / ( vals.max() - vals.min() )

def createFrame( index ):

    print( '  Creating frame {:}/{:}' \
           .format( index+1, FileNumberArray.shape[0] ) )

    ds = yt.load( PlotDirectory + '{:}.plt{:}' \
                  .format( ProblemName, \
                           str( FileNumberArray[index] ).zfill( 8 ) ) )

    im, sc = yt.volume_render( ds, ('boxlib','PF_E') )

    cam = sc.get_camera()
    x = 1.7
    cam.set_width( [ x, x, x ] )
#    cam.set_position( [ 1.0, 0.0, 0.0 ] )
#    cam.set_focus( [ 0.0, 0.0, 0.0 ] )

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

    fileName = '{:}.png'.format( str( index ).zfill( 3 ) )
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

    print( '  Saved {:}'.format( fileName ) )

#createFrame( 0 )
for i in range( 0, FileNumberArray.shape[0], 1 ):
    createFrame( i )

'''
To make a movie, generate the *.png files
(be sure to remove the *Render*png files)
and exectute
ffmpeg -framerate 30 -pattern_type glob -i '*png' -c:v libx264 -pix_fmy yuv420p out.mp4
'''
