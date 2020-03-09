#!/usr/local/bin/python3

import yt
import numpy as np
import subprocess
from os import listdir
from os.path import isfile
from sys import argv, exit
import matplotlib.pyplot as plt

"""
Default use (python3 PlotAMReX_yt.py), plots last plot-file in ProblemDirectory

Can also specify certain plot-file:  python3 PlotAMReX_yt.py thornado_00000010
"""

yt.funcs.mylog.setLevel(0) # Suppress initial yt output to screen

# --- Get user's HOME directory ---
HOME = subprocess.check_output( ["echo $HOME"], shell = True)
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# --- Get user's THORNADO_DIR directory ---
THORNADO_DIR = subprocess.check_output( ["echo $THORNADO_DIR"], shell = True)
THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

#### ========== User Input ==========

# Specify directory containing plotfiles
ProblemDirectory = THORNADO_DIR + 'SandBox/AMReX/TwoMoment_Test/'

# Specify name of problem (only used for name of output file(s))
ProblemName = 'KHI'

# Specify plot file base name
PlotFileBaseName = 'thornado'

# Specify field to plot
VariableToPlot = 'CR_N'

# Specify to plot in log-scale
UseLogScale = False

# Specify whether or not to use physical units
UsePhysicalUnits = False

# Specify coordinate system (currently supports 'cartesian' and 'spherical' )
CoordinateSystem = 'cartesian'

# Specify aspect ratio (relativistic KHI needs aspect = 0.5)
aspect = 1.0

# Specify colormap
cmap = 'jet'

# Specify to write out data to a file
WriteOut = 1
# Specify to read in data from a file
ReadIn = 1

#init or final
time = 1
#number of cells
cells = "8"

# Specify whether or not to make a movie
MakeMovie, DataFileName, TimeFileName = False, 'MovieData.dat', 'MovieTime.dat'

#### ================================

if( len( argv ) > 1 ):
    File = argv[1]
else:
    # Get last plotfile in directory
    FileArray \
      = np.sort(np.array( [ file for file in listdir( ProblemDirectory ) ] ))
    FileList = []
    for iFile in range( FileArray.shape[0] ):
        sFile = FileArray[iFile]
        if( sFile[0:len(PlotFileBaseName)+1] == PlotFileBaseName + '_' \
              and sFile[len(PlotFileBaseName)+1].isdigit() ):
            FileList.append( sFile )
    FileArray = np.array( FileList )
    File = FileArray[-1]

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
if  ( VariableToPlot == 'PR_D'  ):
    Data = CoveringGrid['PR_D_001_001' ].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'g/cm**3'
elif( VariableToPlot == 'PR_I1' ):
    Data = CoveringGrid['PR_I1_001_001'].to_ndarray()
    if( UsePhysicalUnits ):
        if  ( CoordinateSystem == 'cartesian' ):
            DataUnit = 'km/s'
        elif( CoordinateSystem == 'spherical' ):
            DataUnit = 'km/s'
elif( VariableToPlot == 'PR_I2' ):
    Data = CoveringGrid['PR_I2_001_001'].to_ndarray()
    if( UsePhysicalUnits ):
        if  ( CoordinateSystem == 'cartesian' ):
            DataUnit = 'km/s'
        elif( CoordinateSystem == 'spherical' ):
            DataUnit = '1/s'
elif( VariableToPlot == 'PR_I3' ):
    Data = CoveringGrid['PR_I3_001_001'].to_ndarray()
    if( UsePhysicalUnits ):
        if  ( CoordinateSystem == 'cartesian' ):
            DataUnit = 'km/s'
        elif( CoordinateSystem == 'spherical' ):
            DataUnit = '1/s'
elif( VariableToPlot == 'CR_N'  ):
    Data = CoveringGrid['CR_N _001_001' ].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'g/cm**3'
elif( VariableToPlot == 'CR_G1' ):
    Data = CoveringGrid['CR_G1_001_001'].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'g/cm**2/s'
elif( VariableToPlot == 'CR_G2' ):
    Data = CoveringGrid['CR_G2_001_001'].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'g/cm**2/s'
elif( VariableToPlot == 'CR_I3' ):
    Data = CoveringGrid['CR_I3_001_001'].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'g/cm**2/s'

    # Apply boundary conditions to theta elements
    for j in indX1:

        # North pole
        k = 0
        Data[j,k,0] \
          = 1.0 / X1[j] \
              * ( ( X1[j+1]**2 * PF_V2[j+1,k] \
                      - X1[j-1]**2 * PF_V2[j-1,k] ) \
                  / ( 2.0 * dX1 ) \
                - ( PF_V1[j,k+1] \
                      - PF_V1[j,k] ) \
                  / ( 2.0 * dX2 ) )

        # South pole
        k = nX[1]-1
        Data[j,k,0] \
          = 1.0 / X1[j] \
              * ( ( X1[j+1]**2 * PF_V2[j+1,k] \
                      - X1[j-1]**2 * PF_V2[j-1,k] ) \
                  / ( 2.0 * dX1 ) \
                - ( PF_V1[j,k] \
                      - PF_V1[j,k-1] ) \
                  / ( 2.0 * dX2 ) )

    if( UsePhysicalUnits ):
        DataUnit = '1/s'


if  (WriteOut == 1 ):
    if (time == 1):
        f = open("thornado" + cells + "fin" + ".txt", "w")
        f.write("\n".join([",".join([str(n) for n in item]) for item in Data[:,:,0].tolist()]))
        f.close()
    if (time == 0):
        f = open("thornado" + cells + "init" + ".txt", "w")
        f.write("\n".join([",".join([str(n) for n in item]) for item in Data[:,:,0].tolist()]))
        f.close()

#if  ( ReadIn == 1 ):
#    dataList = [[float(s) for s in item.strip().split(",")] for item in open("numpytest.txt").readlines()]
#    rows = len(dataList)
#    cols = len(dataList[0])
#    arr1 = np.zeros([rows,cols], dtype=np.float64)
#    for i, row in enumerate(dataList):
#        for j, number in enumerate(row):
#            arr1[i][j] = number
#

if  ( nDims == 1 ):

    x = np.linspace( xL[0].to_ndarray(), xH[0].to_ndarray(), nX[0] )

    plt.plot( x, Data[:,0,0], 'k-' )
    plt.xlim( xL[0], xH[0] )
    plt.xlabel( 'X1' )
    plt.ylabel( VariableToPlot )
    plt.show()
    plt.close()

    if( MakeMovie ):
        print( 'Making a movie...' )
        from matplotlib import animation

        Overwrite = True
        if( isfile( DataFileName ) ):
            Overwrite = input( 'File: "{:}" exists. overwrite? (Y/N): '.format \
                          ( DataFileName ) )
            if( not Overwrite == 'Y' ):
                print( 'Not overwriting file, using existing file for movie.' )
                Overwrite = False
            else:
                Overwrite = True

        if( Overwrite ):
            # Put all time-slices into one array to use for movie making
            Data = np.empty( (FileArray.shape[0]+1,nX[0]), float )
            Time = np.empty( FileArray.shape[0], float )

            Data[0] = x
            for i in range( FileArray.shape[0] ):
                print( '{:}/{:}'.format( i+1, FileArray.shape[0] ) )
                ds = yt.load( '{:}'.format( ProblemDirectory + FileArray[i] ) )

                CoveringGrid \
                  = ds.covering_grid \
                      ( level           = MaxLevel, \
                        left_edge       = xL, \
                        dims            = nX * 2**MaxLevel, \
                        num_ghost_zones = nX[0] )

                Data[i+1] = CoveringGrid[VariableToPlot].to_ndarray()[:,0,0]
                Time[i] = ds.current_time

            np.savetxt( DataFileName, Data )
            np.savetxt( TimeFileName, Time )

        Data = np.loadtxt( DataFileName, skiprows=1 )
        Time = np.loadtxt( TimeFileName )

        fig, ax = plt.subplots()

        ax.set_xlim( xL[0], xH[0] )
        ax.set_ylim( np.min( Data ), np.max( Data ) )

        ax.set_xlabel( 'X1' )
        ax.set_ylabel( VariableToPlot )

        Width     = xH[0] - xL[0]
        Height    = np.max( Data ) - np.min( Data )
        time_text = plt.text( xL[0] + 0.5 * Width, \
                              np.min( Data ) + 0.7 * Height, '' )

        if( UseLogScale ): ax.set_yscale( 'log' )
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
    Data = np.loadtxt( DataFileName )
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
                        aspect = aspect, \
                        origin = 'lower-left-window' )

    slc.set_cmap( field = field, cmap = cmap )

    slc.set_log( field, UseLogScale )
    #slc.set_zlim( field, 0.0, 2.0 ) # Set colorbar limits
    #slc.set_colorbar_label( field, 'Primitive Rest-Mass-Density' )

    if( CoordinateSystem == 'spherical' ):
        slc.set_width( 2 * xH[0].to_ndarray(), length_unit )

    slc.save( ProblemName + '_' + PlotFileBaseName + '_' + VariableToPlot \
                + '_{:}.png'.format( File[-8:] ) )

    if( MakeMovie ):
        from matplotlib import animation

        Overwrite = True
        if( isfile( DataFileName ) ):
            Overwrite = input( 'File: "{:}" exists. overwrite? (Y/N): '.format \
                          ( DataFileName ) )
            if( not Overwrite == 'Y' ):
                print( 'Not overwriting file, using existing file for movie.' )
                Overwrite = False
            else:
                Overwrite = True

        if( Overwrite ):
            # Put all time-slices into one array to use for movie making
            Data = np.empty( (FileArray.shape[0],nX[0],nX[1]), float )
            Time = np.empty( FileArray.shape[0], float )
            for i in range( FileArray.shape[0] ):
                print( '{:}/{:}'.format( i+1, FileArray.shape[0] ) )
                ds = yt.load( '{:}'.format( ProblemDirectory + FileArray[i] ) )

                CoveringGrid \
                  = ds.covering_grid \
                      ( level           = MaxLevel, \
                        left_edge       = xL, \
                        dims            = nX * 2**MaxLevel, \
                        num_ghost_zones = nX[0] )

                Data[i] = CoveringGrid[VariableToPlot].to_ndarray()[:,:,0]
                Time[i] = ds.current_time

            # Save multi-D array with np.savetxt. Taken from:
            # https://stackoverflow.com/questions/3685265/how-to-write-a-multidimensional-array-to-a-text-file

            with open( DataFileName, 'w' ) as FileOut:
                FileOut.write( '# Array shape: {0}\n'.format( Data.shape ) )

                # Iterating through an n-dimensional array produces slices along
                # the last axis. This is equivalent to Data[i] in this case
                for TimeSlice in Data:
                    np.savetxt( FileOut, TimeSlice )
                    FileOut.write( '# New slice\n' )

            np.savetxt( 'MovieTime.dat', Time )

        Data = np.loadtxt( DataFileName ).reshape( \
                                            (FileArray.shape[0],nX[0],nX[1]) )
        Time = np.loadtxt( 'MovieTime.dat' )

        fig = plt.figure()
        def f(t):
            return Data[t].T

        PlotTranspose = False
        if( PlotTranspose ):
            fig.suptitle( '|{:}-{:}.T|'.format \
                          ( VariableToPlot, VariableToPlot ), \
                          fontsize = 20 )
            def f(t):
                return np.abs( Data[t] - Data[t].T )

        if( UseLogScale ):
            from matplotlib.colors import LogNorm
            vmin = min( +np.inf, np.min( np.abs( Data ) ) )
            vmax = max( -np.inf, np.max( np.abs( Data ) ) )

            im = plt.imshow( np.abs(f(0)), cmap = cmap, animated = True, \
                             vmin = vmin, vmax = vmax, \
                             extent = [ xL[0], xH[0], xL[1], xH[1] ], \
                             aspect = 'equal', \
                             origin = 'lower', norm = LogNorm() )
        else:
            im = plt.imshow( f(0), \
                             cmap = cmap, \
                             animated = True, \
                             vmin = np.min(Data), vmax = np.max(Data), \
                             extent = [ xL[0], xH[0], xL[1], xH[1] ], \
                             aspect = 'equal', \
                             origin = 'lower'  )

        plt.colorbar(im)

        Width  = xH[0] - xL[0]
        Height = xH[1] - xL[1]

        time_text = plt.text( xL[0] + 0.5 * Width, xL[1] + 0.9 * Height, '' )

        if( UseLogScale ):
            def UpdateFrame(t):
                im.set_array( np.abs(f(t)) )
                time_text.set_text('time = {:.3e}'.format( Time[t] ) )
                return im,
        else:
            def UpdateFrame(t):
                im.set_array( f(t) )
                time_text.set_text('time = {:.3e}'.format( Time[t] ) )
                return im,

        # Call the animator
        anim = animation.FuncAnimation \
                 ( fig, UpdateFrame, frames = FileArray.shape[0], \
                   interval = 100, blit = True)

        anim.save( '{:}_{:}.mp4'.format( ProblemName, VariableToPlot ), \
                   fps = 5 )
        plt.close()

    exit()

elif ( nDims == 3 ):
    print( 'nDims = 3 not yet implemented. Exiting...\n' )
    exit()
else:
    print( 'Invalid choice of nDims ({:}). Exiting...\n'.format( nDims ) )
    exit()
