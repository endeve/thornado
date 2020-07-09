#!/usr/local/bin/python3

import yt
import numpy as np
import subprocess
from os import listdir
from os.path import isfile
from sys import argv, exit
import matplotlib.pyplot as plt

"""
Default usage, plots last plotfile in DataDirectory:

  $ python3 PlotFieldsAMReX.py

Alernate usage, plot specific file in DataDirectory:

  $ python3 PlotFieldsAMReX.py thornado_00000010

TODO:
  - Add SymLogNorm scaling
"""

# https://yt-project.org/doc/faq/index.html#how-can-i-change-yt-s-log-level
yt.funcs.mylog.setLevel(40) # Suppress yt warnings

# --- Get user's HOME directory ---

HOME = subprocess.check_output( ["echo $HOME"], shell = True )
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# --- Get user's THORNADO_DIR directory ---

THORNADO_DIR = subprocess.check_output( ["echo $THORNADO_DIR"], shell = True )
THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

#### ========== User Input ==========

# Specify name of problem
ProblemName = 'SASI'

# Specify directory containing plotfiles
DataDirectory = THORNADO_DIR + 'SandBox/AMReX/'

# Specify plot file base name
PlotFileBaseName = 'thornado'

# Specify field to plot
Field = 'PF_D'

# Specify to plot in log-scale
UseLogScale = True

# Specify whether or not to use physical units
UsePhysicalUnits = True

# Specify coordinate system (currently supports 'cartesian' and 'spherical' )
CoordinateSystem = 'spherical'

# Specify aspect ratio (relativistic KHI needs aspect = 0.5)
aspect = 1.0

# Specify colormap
cmap = 'Purples'

# Specify whether or not to make a datafile and/or movie
MakeDataFile = False
MakeMovie    = False
ID = '{:}_{:}_{:}'.format( ProblemName, PlotFileBaseName[4:], Field )
DataFileName = '{:}_MovieData.dat'.format( ID )
TimeFileName = '{:}_MovieTime.dat'.format( ID )

#### ====== End of User Input =======

# Append "/" to DataDirectory, if not present
if( not DataDirectory[-1] == '/' ): DataDirectory += '/'

if( len( argv ) == 1 ):

    # Get last plotfile in directory

    FileArray \
      = np.sort( np.array( [ file for file in listdir( DataDirectory ) ] ) )

    FileList = []

    for iFile in range( FileArray.shape[0] ):

        sFile = FileArray[iFile]

        if( sFile[0:len(PlotFileBaseName)+1] == PlotFileBaseName + '_' \
              and sFile[len(PlotFileBaseName)+1].isdigit() ):
            FileList.append( sFile )

    FileArray = np.array( FileList )
    File      = FileArray[-1]

elif( len( argv ) == 2 ):

    File = argv[1]

else:

    print( 'Invalid number of optional parameters' )
    exit( 'Exiting...' )

# Remove "/" at end of filename, if present
if ( File[-1] == '/' ): File = File[:-1]

ds = yt.load( '{:}'.format( DataDirectory + File ) )

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

TimeUnit = ''
if( UsePhysicalUnits ): TimeUnit = 'ms'

"""
https://yt-project.org/doc/reference/api/
yt.data_objects.construction_data_containers.html#yt.data_objects.
construction_data_containers.YTCoveringGrid
"""
CoveringGrid \
  = ds.covering_grid \
      ( level           = MaxLevel, \
        left_edge       = xL, \
        dims            = nX * 2**MaxLevel, \
        num_ghost_zones = nX[0] )

# XXX.to_ndarray() strips yt array of units

DataUnit = ''

if  ( Field == 'PF_D'  ):
    Data = CoveringGrid['PF_D' ].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'g/cm**3'

elif( Field == 'PF_V1' ):
    Data = CoveringGrid['PF_V1'].to_ndarray()
    if( UsePhysicalUnits ):
        if  ( CoordinateSystem == 'cartesian' ):
            DataUnit = 'km/s'
        elif( CoordinateSystem == 'spherical' ):
            DataUnit = 'km/s'

elif( Field == 'PF_V2' ):
    Data = CoveringGrid['PF_V2'].to_ndarray()
    if( UsePhysicalUnits ):
        if  ( CoordinateSystem == 'cartesian' ):
            DataUnit = 'km/s'
        elif( CoordinateSystem == 'spherical' ):
            DataUnit = '1/s'

elif( Field == 'PF_V3' ):
    Data = CoveringGrid['PF_V3'].to_ndarray()
    if( UsePhysicalUnits ):
        if  ( CoordinateSystem == 'cartesian' ):
            DataUnit = 'km/s'
        elif( CoordinateSystem == 'spherical' ):
            DataUnit = '1/s'

elif( Field == 'PF_E'  ):
    Data = CoveringGrid['PF_E' ].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'erg/cm**3'

elif( Field == 'CF_D'  ):
    Data = CoveringGrid['CF_D' ].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'g/cm**3'

elif( Field == 'CF_S1' ):
    Data = CoveringGrid['CF_S1'].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'g/cm**2/s'

elif( Field == 'CF_S2' ):
    Data = CoveringGrid['CF_S2'].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'g/cm**2/s'

elif( Field == 'CF_S3' ):
    Data = CoveringGrid['CF_S3'].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'g/cm**2/s'

elif( Field == 'CF_E'  ):
    Data = CoveringGrid['CF_E' ].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'erg/cm**3'

elif( Field == 'AF_P'  ):
    Data = CoveringGrid['AF_P' ].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'erg/cm**3'

elif( Field == 'AF_Cs' ):
    Data = CoveringGrid['AF_Cs'].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'km/s'

elif( Field == 'DF_Sh_X1'  ):
    Data = CoveringGrid['DF_Sh_X1'].to_ndarray()

elif( Field == 'DF_Sh_X2'  ):
    Data = CoveringGrid['DF_Sh_X2'].to_ndarray()

elif( Field == 'DF_Sh_X3'  ):
    Data = CoveringGrid['DF_Sh_X3'].to_ndarray()

# --- Derived Fields ---

elif( Field == 'PolytropicConstant' ):
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

elif( Field == 'Vorticity' ):

    dX1 = ( xH[0].to_ndarray() - xL[0].to_ndarray() ) / nX[0]
    dX2 = ( xH[1].to_ndarray() - xL[1].to_ndarray() ) / nX[1]

    XL = xL.to_ndarray() + 0.5 * np.array( [ dX1, dX2, 0.0 ] )
    XH = xH.to_ndarray() - 0.5 * np.array( [ dX1, dX2, 0.0 ] )

    X1 = np.linspace( XL[0], XH[0], nX[0] )
    X2 = np.linspace( XL[1], XH[1], nX[1] )

    PF_V1 = CoveringGrid['PF_V1'].to_ndarray()
    PF_V2 = CoveringGrid['PF_V2'].to_ndarray()
    indX1 = np.linspace( 1, nX[0]-2, nX[0]-2, dtype = int )
    indX2 = np.linspace( 1, nX[1]-2, nX[1]-2, dtype = int )
    Data = np.zeros( (nX[0],nX[1],1), float )
    for j in indX1:
        for k in indX2:
            Data[j,k,0] \
              = 1.0 / X1[j] \
                  * ( ( X1[j+1]**2 * PF_V2[j+1,k] \
                          - X1[j-1]**2 * PF_V2[j-1,k] ) \
                      / ( 2.0 * dX1 ) \
                    - ( PF_V1[j,k+1] \
                          - PF_V1[j,k-1] ) \
                      / ( 2.0 * dX2 ) )

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
else:

    print( 'Invalid field: {:}'.format( Field ) )
    exit( 'Exiting...' )

if  ( nDims == 1 ):

    dX1 = ( xH[0].to_ndarray() - xL[0].to_ndarray() ) / nX[0]

    x = np.linspace( xL[0].to_ndarray() + 0.5 * dX1, \
                     xH[0].to_ndarray() - 0.5 * dX1, nX[0] )

    plt.plot( x, Data[:,0,0], 'k-' )
    if( UseLogScale ): plt.yscale( 'log' )
    plt.xlim( xL[0], xH[0] )
    plt.xlabel( 'X1' )
    plt.ylabel( Field )
    plt.show()
    plt.close()

    if( MakeDataFile ):

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
                ds = yt.load( '{:}'.format( DataDirectory + FileArray[i] ) )

                CoveringGrid \
                  = ds.covering_grid \
                      ( level           = MaxLevel, \
                        left_edge       = xL, \
                        dims            = nX * 2**MaxLevel, \
                        num_ghost_zones = nX[0] )

                if( Field == 'PolytropicConstant' ):

                    AF_P  = CoveringGrid['AF_P' ].to_ndarray()[:,0,0]
                    AF_Gm = CoveringGrid['AF_Gm'].to_ndarray()[:,0,0]
                    PF_D  = CoveringGrid['PF_D' ].to_ndarray()[:,0,0]
                    Data[i+1] = AF_P / PF_D**(AF_Gm)

                else:

                    Data[i+1] = CoveringGrid[Field].to_ndarray()[:,0,0]

                Time[i] = ds.current_time

            np.savetxt( DataFileName, Data )
            np.savetxt( TimeFileName, Time )

    if( MakeMovie ):

        from matplotlib import animation

        if( not isfile( DataFileName ) ):

            print( 'File: {:} does not exist.'.format( DataFileName ) )
            exit( 'Exiting...' )

        Data = np.loadtxt( DataFileName, skiprows = 1 )
        Time = np.loadtxt( TimeFileName )

        fig, ax = plt.subplots()

        ax.set_xlim( xL[0], xH[0] )
        ax.set_ylim( np.min( Data ), np.max( Data ) )

        ax.set_xlabel( 'X1' )
        ax.set_ylabel( Field )

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
            y = Data[t]
            line.set_data( x, y )
            time_text.set_text( 'time = {:.3e} {:}'.format \
              ( Time[t], TimeUnit ) )
            return line,

        anim = animation.FuncAnimation( fig, UpdateFrame, \
                                        init_func = InitializeFrame, \
                                        frames = FileArray.shape[0], \
                                        interval = 200, blit = True )
        anim.save( '{:}_Movie.mp4'.format( ID ) )

    """
    # Plot slices from 'Data' file, created with MakeMovie script above
    Data = np.loadtxt( DataFileName )
    x = Data[0]
    Data = np.copy( Data[1:] )
    plt.plot( x, Data[0],  'k-',  label = 't = 0 ms'   )
    plt.xlabel( 'Radial Distance [km]' )
    if( UsePhysicalUnits ):
        plt.ylabel( '{:} [{:}]'.format( Field, DataUnit ) )
    else:
        plt.ylabel( '{:}'.format( Field ) )
    plt.xlim( x.min(), x.max() )
    plt.ylim( Data.min(), Data.max() )
    plt.legend()
    plt.show()
    plt.close()
    """

elif( nDims == 2 ):

    '''
    # To make lineout plot
    # From: https://yt-project.org/doc/visualizing/
    #       manual_plotting.html#line-plots

    x = np.linspace( xL[0], xH[0], nX[0] )

    oray = ds.ortho_ray( axis = 0, coords = (0,0) )
    srt  = np.argsort( oray[Field] )

    plt.plot( np.array( oray[Field][srt[::-1]] ) )

    if( UseLogScale ): plt.yscale( 'log' )
    plt.show()
    exit()
    '''

    data          = { Field: (Data,DataUnit) }
    field         = Field
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
                      [ [xL[0],xH[0]], \
                        [xL[1],xH[1]], \
                        [xL[2],xH[2]] ] ), \
             length_unit = length_unit, \
             geometry = CoordinateSystem )

    if( CoordinateSystem == 'spherical' ):

        w0 = xH[0] + xL[0]
        w1 = 2.0 * xH[0]
        slc = yt.SlicePlot( ds        = ds, \
                            normal    = SliceVariable, \
                            fields    = field, \
                            width     = ((w0),(w1)), \
                            axes_unit = length_unit, \
                            aspect    = aspect, \
                            origin    = 'lower-left-window' )

    else:

        slc = yt.SlicePlot( ds        = ds, \
                            normal    = SliceVariable, \
                            fields    = field, \
                            axes_unit = length_unit, \
                            aspect    = aspect, \
                            origin    = 'lower-left-window' )

    slc.set_cmap( field = field, cmap = cmap )

    slc.set_log( field, UseLogScale )
    #slc.set_zlim( field, 0.0, 2.0 ) # Set colorbar limits
    #slc.set_colorbar_label( field, 'Primitive Rest-Mass-Density' )

    slc.save( ProblemName + '_' + PlotFileBaseName + '_' + Field \
                + '_{:}.png'.format( File[-8:] ) )

    if( MakeDataFile ):

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
                ds = yt.load( '{:}'.format( DataDirectory + FileArray[i] ) )

                CoveringGrid \
                  = ds.covering_grid \
                      ( level           = MaxLevel, \
                        left_edge       = xL, \
                        dims            = nX * 2**MaxLevel, \
                        num_ghost_zones = nX[0] )

                Data[i] = CoveringGrid[Field].to_ndarray()[:,:,0]
                Time[i] = ds.current_time

            """
            Save multi-D array with np.savetxt. Taken from:
            https://stackoverflow.com/questions/3685265/
            how-to-write-a-multidimensional-array-to-a-text-file
            """
            with open( DataFileName, 'w' ) as FileOut:
                FileOut.write( '# Array shape: {0}\n'.format( Data.shape ) )

                # Iterating through an n-dimensional array produces slices along
                # the last axis. This is equivalent to Data[i] in this case
                for TimeSlice in Data:
                    np.savetxt( FileOut, TimeSlice )
                    FileOut.write( '# New slice\n' )

            np.savetxt( TimeFileName, Time )

        if( MakeMovie ):

            from matplotlib import animation

            if( not isfile( DataFileName ) ):

                print( 'File: {:} does not exist.'.format( DataFileName ) )
                exit( 'Exiting...' )

            Data \
              = np.loadtxt( DataFileName ).reshape( \
                  (FileArray.shape[0],nX[0],nX[1]) )
            Time = np.loadtxt( TimeFileName )

            fig = plt.figure()
            def f(t):
                return Data[t].T

            PlotTranspose = False
            if( PlotTranspose ):
                fig.suptitle( '|{:}-{:}.T|'.format \
                              ( Field, Field ), \
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

            time_text \
              = plt.text( xL[0] + 0.5 * Width, xL[1] + 0.9 * Height, '' )

            if( UseLogScale ):
                def UpdateFrame(t):
                    im.set_array( np.abs(f(t)) )
                    time_text.set_text( 'Time = {:.3e} {:}'.format \
                      ( Time[t], TimeUnit ) )
                    return im,
            else:
                def UpdateFrame(t):
                    im.set_array( f(t) )
                    time_text.set_text( 'Time = {:.3e} {:}'.format \
                      ( Time[t], TimeUnit ) )
                    return im,

            # Call the animator
            anim = animation.FuncAnimation \
                     ( fig, UpdateFrame, frames = FileArray.shape[0], \
                       interval = 100, blit = True)

            anim.save( '{:}_Movie.mp4'.format( ID ), fps = 5 )
            plt.close()

    exit()

elif ( nDims == 3 ):
    print( 'nDims = 3 not yet implemented. Exiting...\n' )
    exit()
else:
    print( 'Invalid choice of nDims ({:}). Exiting...\n'.format( nDims ) )
    exit()
