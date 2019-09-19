#!/usr/local/bin/python3

import yt
import numpy as np
import subprocess
from os import listdir
from sys import argv, exit
import matplotlib.pyplot as plt

"""
To-Do:
  - Make PNS circular
  - Add 2D movie-making capability
"""

yt.funcs.mylog.setLevel(0) # Suppress initial yt output to screen

# --- Get user's HOME directory ---
HOME = subprocess.check_output( ["echo $HOME"], shell = True)
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# Specify directory containing plotfiles

#'''
ProblemDirectory \
  = HOME + 'Research/SN/thornado/SandBox/AMReX/Euler_Relativistic/'
ProblemName = 'SAS_R'
VariableToPlot = 'Entropy'
#'''

'''
ProblemDirectory \
  = HOME + 'Desktop/'
ProblemName = 'SAS_R'
VariableToPlot = 'PF_D'
'''

UsePhysicalUnits = False
CoordinateSystem = 'cartesian'
if( ProblemName == 'SAS' ):
    CoordinateSystem = 'spherical'
if( ProblemName == 'SAS_R' ):
    CoordinateSystem = 'spherical'
    UsePhysicalUnits = True

# Get last plotfile in directory
FileArray \
  = np.sort(np.array( [ file for file in listdir( ProblemDirectory ) ] ))
FileList = []
for iFile in range( FileArray.shape[0] ):
    sFile = FileArray[iFile]
    if( sFile[0:8] == 'thornado' ):
        FileList.append( sFile )
FileArray = np.array( FileList )
File = FileArray[-1]

# Also allow user to plot specific file
if( len( argv ) > 1 ):
    File = argv[1]

# Change aspect ratio to make plotting cells square
aspect = 1.0
if( ProblemName == 'KHI_R' ): aspect = 0.5

# Remove "/" at end of filename, if present
if ( File[-1] == '/' ): File = File[:-1]

ds = yt.load( '{:}'.format( ProblemDirectory + File ) )

print( 'Reading from file: {:}'.format( File ) )
MaxLevel = ds.index.max_level
Time     = ds.current_time
nX       = ds.domain_dimensions
xL       = ds.domain_left_edge
xH       = ds.domain_right_edge

# Get dimensionality of problem
if  ( nX[1] == 1 and nX[2] == 1 ):
    nDims = 1
elif( nX[2] == 1 ):
    nDims = 2
else:
    nDims = 3

CoveringGrid \
  = ds.covering_grid \
      ( level           = MaxLevel, \
        left_edge       = xL, \
        dims            = nX * 2**MaxLevel, \
        num_ghost_zones = nX[0] )

# XXX.to_ndarray() strips array of yt units

DataUnit = ''
if  ( VariableToPlot == 'PF_D'  ):
    Data = CoveringGrid['PF_D' ].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'g/cm**3'
elif( VariableToPlot == 'PF_V1' ):
    Data = CoveringGrid['PF_V1'].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'km/s'
elif( VariableToPlot == 'PF_V2' ):
    Data = CoveringGrid['PF_V2'].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'km/s'
elif( VariableToPlot == 'PF_V3' ):
    Data = CoveringGrid['PF_V3'].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'km/s'
elif( VariableToPlot == 'PF_E'  ):
    Data = CoveringGrid['PF_E' ].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'erg/cm**3'
elif( VariableToPlot == 'CF_D'  ):
    Data = CoveringGrid['CF_D' ].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'g/cm**3'
elif( VariableToPlot == 'CF_S1' ):
    Data = CoveringGrid['CF_S1'].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'g/cm**2/s'
elif( VariableToPlot == 'CF_S2' ):
    Data = CoveringGrid['CF_S2'].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'g/cm**2/s'
elif( VariableToPlot == 'CF_S3' ):
    Data = CoveringGrid['CF_S3'].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'g/cm**2/s'
elif( VariableToPlot == 'CF_E'  ):
    Data = CoveringGrid['CF_E' ].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'erg/cm**3'
elif( VariableToPlot == 'AF_P'  ):
    Data = CoveringGrid['AF_P' ].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'erg/cm**3'
elif( VariableToPlot == 'Entropy' ):
    PF_D  = CoveringGrid['PF_D' ].to_ndarray()
    AF_P  = CoveringGrid['AF_P' ].to_ndarray()
    AF_Gm = CoveringGrid['AF_Gm'].to_ndarray()
    Data  = AF_P / PF_D**AF_Gm
    if( UsePhysicalUnits ):
        if( round( AF_Gm[0][0][0], 1 ) == 1.4 ):
            DataUnit = 'erg/cm**3/(g/cm**3)**({:})'.format( \
                         AF_Gm[0][0][0] )
        else:
            DataUnit = 'erg/cm**3/(g/cm**3)**({:}/3)'.format( \
                         int( 3 * AF_Gm[0][0][0] ) )

if  ( nDims == 1 ):

    x = np.linspace( xL[0].to_ndarray(), xH[0].to_ndarray(), nX[0] )

    #"""
    plt.plot( x, Data[:,0,0], 'k-' )
    plt.xlim( xL[0], xH[0] )
    plt.xlabel( 'X1' )
    plt.ylabel( VariableToPlot )
    plt.show()
    plt.close()
    #"""

    MakeMovie = False
    if( MakeMovie ):
        print( 'Making a movie...' )
        from matplotlib import animation

        #'''
        # Put all time-slices into one array to use for movie making
        Data = np.empty( (FileArray.shape[0],nX[0]), float )
        Time = np.empty( FileArray.shape[0], float )
        for i in range( FileArray.shape[0] ):
            print( '{:}/{:}'.format( i, FileArray.shape[0] ) )
            ds = yt.load( '{:}'.format( ProblemDirectory + FileArray[i] ) )

            CoveringGrid \
              = ds.covering_grid \
                  ( level           = MaxLevel, \
                    left_edge       = xL, \
                    dims            = nX * 2**MaxLevel, \
                    num_ghost_zones = nX[0] )

            Data[i] = CoveringGrid[VariableToPlot].to_ndarray()[:,0,0]
            Time[i] = ds.current_time

        np.savetxt( 'Data', Data )
        np.savetxt( 'Time', Time )
        #'''
        Data = np.loadtxt( 'Data' )
        Time = np.loadtxt( 'Time' )

        fig, ax = plt.subplots()

        ax.set_xlim( xL[0], xH[0] )
        ax.set_ylim( np.min( Data ), np.max( Data ) )

        ax.set_xlabel( 'X1' )
        ax.set_ylabel( VariableToPlot )

        Width     = xH[0] - xL[0]
        Height    = np.max( Data ) - np.min( Data )
        time_text = plt.text( xL[0] + 0.5 * Width, \
                              np.min( Data ) + 0.7 * Height, '' )

        ax.set_yscale( 'log' )
        #IC,   = ax.plot([],[], color = 'red',   linewidth = 2 )
        line, = ax.plot([],[], color = 'black', linewidth = 1 )
        def InitializeFrame():
            #IC.set_data([],[])
            line.set_data([],[])
            time_text.set_text('')
            return line, time_text#, IC
        def UpdateFrame(t):
            #IC.set_data( x, Data[0] )
            y = np.abs( Data[t] )
            line.set_data( x, y )
            if( UsePhysicalUnits ):
                time_text.set_text( 'time = {:d} ms'.format( int( Time[t] ) ) )
            else:
                time_text.set_text( 'time = {:d}'.format( int( Time[t] ) ) )
            return line,

        anim = animation.FuncAnimation( fig, UpdateFrame, \
                                        init_func = InitializeFrame, \
                                        frames = FileArray.shape[0], \
                                        interval = 200, blit = True )
        anim.save( VariableToPlot + '.mp4' )

    """
    # Plot slices from 'Data' file, created with MakeMovie script above
    Data = np.loadtxt( 'Data' )
    plt.semilogy( x, Data[0],  'k-',  label = 't = 0 ms'   )
    plt.semilogy( x, Data[-1], 'k--', label = 't = 300 ms' )
    plt.xlabel( 'Radial Distance [km]' )
    plt.ylabel( r'$\rho\ \left[g\,cm^{-3}\right]$' )
    #plt.suptitle( 'Standing Accretion Shock' )
    plt.xlim( 40, 360 )
    plt.ylim( 8.0e6, 3.0e10 )
    plt.text( 160, 1e10, 'nX = 32', fontsize = 20 )
    plt.legend()
    plt.show()
    plt.close()
    """

elif( nDims == 2 ):

    '''
    # To make lineout plot
    # From: https://yt-project.org/doc/reference/api/yt.data_objects.selection_data_containers.html#yt.data_objects.selection_data_containers.YTOrthoRay
    oray = ds.ortho_ray( axis = 0, coords = (0,0) )
    x = np.linspace( xL[0], xH[0], nX[0] )
    plt.plot( x, oray[VariableToPlot] )
    plt.show()
    exit()
    '''

    data          = { VariableToPlot: (Data,DataUnit) }
    field         = VariableToPlot
    length_unit   = 'code_length'
    SliceVariable = 'z'
    if( CoordinateSystem == 'spherical' ):
        SliceVariable = 'phi'
    if( UsePhysicalUnits ):
        length_unit = 'km'

    ds = yt.load_uniform_grid \
           ( data, \
             nX, \
             bbox = np.array( \
                      [ [xL[0],xH[0]], [xL[1],xH[1]], [xL[2],xH[2]] ] ), \
             length_unit = length_unit, \
             geometry = CoordinateSystem )

    slc = yt.SlicePlot( ds, SliceVariable, field, \
                        axes_unit = length_unit, \
                        aspect = aspect )

    slc.set_cmap( field = field, cmap = 'jet' )

    #slc.set_log( field, False ) # set_log is True by default
    #slc.set_zlim( field, 0.0, 2.0 ) # Set colorbar limits
    #slc.set_colorbar_label( field, 'Primitive Rest-Mass-Density' )

    # Only works in Cartesian coordinates
    #slc.annotate_text( [-0.48,0.9], 't = 3.00', coord_system = 'plot', \
    #                   text_args={'color':'black'} )

    if( CoordinateSystem == 'spherical' ):
        slc.set_width( 2 * xH[0].to_ndarray(), length_unit )

    slc.save( ProblemName + '_' + VariableToPlot \
                + '_{:}.png'.format( File[-8:] ) )

    exit()

elif ( nDims == 3 ):
    print( 'nDims = 3 not yet implemented. Exiting...\n' )
    exit()
else:
    print( 'Invalid choice of nDims ({:}). Exiting...\n'.format( nDims ) )
    exit()
