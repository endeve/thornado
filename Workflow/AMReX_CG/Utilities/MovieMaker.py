#!/usr/bin/env python3

import numpy as np
from sys import argv
from matplotlib import animation
from functools import partial
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

import GlobalVariables.Settings as gvS
import GlobalVariables.Units    as gvU
from Utilities.RefinementBoundaryFinder import FindRefinementBoundaries

#=============================================#
#   Included Routines
#
#   MakeMovie
#   CreateMovie
#   CreateFrame         - Contains Line Settings -
#   InitializeFrame
#   UpdateFrame
#   fetchData_AMReX
#   ReadHeader
#   ApplyAction
#
#=============================================#

 #=============================================#
#                                               #
#   MakeMovie                                   #
#                                               #
 #=============================================#
def MakeMovie(  FileNumberArray,    \
                Field,              \
                DataDirectories,    \
                Action = 'None'     ):

    global nDirs
    global nFiles

    # Check if FileNumberArray and DataDirectories are lists
    if type(FileNumberArray) is list:
        nFiles = len(FileNumberArray)
    else:
        nFiles = 1
        
        
    if type(DataDirectories) is list:
        nDirs = len(DataDirectories)
    else:
        nDirs = 1
    
            
    if type(Field) is list:
        nFields = len(Field)
        if nFields == 1:
            Field = Field[0]
    else:
        nFields = 1
    
    
    if  nFiles == nDirs:

        CreateMovie(FileNumberArray,    \
                    Field,              \
                    DataDirectories,    \
                    Action = Action     )
                                
    else:
    
        msg =  "\n MakeMovie Error \n"
        msg += "MakeMovie requires the same number of FileNumberArrays \n"
        msg += "as DataDirectories. One FileNumberArray per DataDirectory \n"
        msg += "must be passed into the routine."
        
        assert(nFiles == nDirs),msg

    return

 #=============================================#
#                                               #
#   CreateMovie                                 #
#                                               #
 #=============================================#
def CreateMovie(FileNumberArray,    \
                Field,              \
                DataDirectory,      \
                Action              ):

    global Lines
    global RefLines
    global time_text
    global IC
    global mesh
    global Data0
    global X1_C0
    global dX10
    global nSS
    global nDirs
    global nFiles
    
    for i in range(nDirs):
        if DataDirectory[i][-1] != '/': DataDirectory[i] += '/'
        
#    Data0, DataUnits, X1_C0, dX10, Time = fetchData_AMReX(-1,                 \
#                                                          FileNumberArray[-1],\
#                                                          DataDirectory[0],  \
#                                                          Field              )

    Data0, DataUnits, X1_C0, dX10, Time = fetchData_AMReX(0,                 \
                                                          FileNumberArray[0],\
                                                          DataDirectory[0],  \
                                                          Field              )
                                                          

    nSS = len(FileNumberArray[0])

    if not gvS.UseCustomLimits:
        gvS.vmin = +np.inf
        gvS.vmax = -np.inf
        for i in range( nFiles ):
            for j in range( nSS ):
                DataFile \
                  = DataDirectory[i] + str(FileNumberArray[i][j]) + '/{:}.dat'.format( Field )
                DataShape, DataUnits, MinVal, MaxVal = ReadHeader( DataFile )
                gvS.vmin = min( gvS.vmin, MinVal )
                gvS.vmax = max( gvS.vmax, MaxVal )
            gvS.vmax = gvS.vmax*1.1

    nX = np.shape(X1_C0)

    xL = X1_C0[0 ] - 0.5 * dX10[0 ]
    xH = X1_C0[-1] + 0.5 * dX10[-1]

    fig = plt.figure()
    ax  = fig.add_subplot( 111 )
    
    
    CreateFrame( ax, xL, xH, dX10, Field, DataUnits )

    anim = animation.FuncAnimation( fig,                                        \
                                    partial(UpdateFrame,                        \
                                            FileNumberArray=FileNumberArray,    \
                                            DataDirectory=DataDirectory,        \
                                            Field = Field,                      \
                                            Action = Action),                   \
                                    init_func = partial(InitializeFrame,        \
                                                        FileNumberArray=FileNumberArray,    \
                                                        DataDirectory=DataDirectory,        \
                                                        Field = Field,                      \
                                                        Action = Action),       \
                                    frames = nSS,                               \
                                    blit = True )



    fps = max( 1, nSS / gvS.MovieRunTime )
    


    print( '\n  Making movie' )
    print( '  ------------' )
    anim.save( gvS.MovieName, fps = fps, dpi = 300 )

    import os
    os.system( 'rm -rf __pycache__ ' )



    return




 #=============================================#
#                                               #
#   CreateFrame                                 #
#                                               #
 #=============================================#
def CreateFrame( ax, xL, xH, dX10, Field, DataUnits ):

    global time_text
    global elem_text
    global Lines
    global RefLines
    global IC
    global mesh


    time_text = ax.text( 0.1, 0.15, '', transform = ax.transAxes, fontsize = 13 )
    elem_text = ax.text( 0.1, 0.1, '', transform = ax.transAxes, fontsize = 13 )

    Lines = ['None']*nDirs
    for i in range(nDirs):
        Lines[i], = ax.plot( [],[],                                         \
                             color  = 'blue',                               \
                             marker = '.',                                  \
                             label  = r'$u_{:}\left(t\right)$'.format(i),   \
                             zorder = 10 )
    
        #color  = plt.cm.Set1(i),                       \
    
    if gvS.PlotMesh:
        mesh,     = ax.plot( [],[] )
    
    
    if gvS.ShowIC:
        IC = ['None']*nDirs
        for i in range(nDirs):
            IC[i], = ax.plot( [],[],                                                \
                              color  = 'k',                                         \
                              linestyle = '--',                                     \
                              label  = r'$u_{{{:},0}}\left(t\right)$'.format(i),    \
                              zorder = 1000,                                        \
                              alpha  = 0.5    )
    
    
    if gvS.ShowRefinement:
        RefLines = ['None']*gvS.RefinementLevels
        for i in range(gvS.RefinementLevels):
            RefLines[i], = ax.plot( [],[],              \
                                    color  = 'red',     \
                                    linestyle = '--',   \
                                    scaley = False,     \
                                    zorder = 0,         \
                                    alpha  = 0.7        )


    ApplyMovieSettings( ax, xL, xH, dX10, Field, DataUnits )

    return






 #=============================================#
#                                               #
#   InitializeFrame                             #
#                                               #
 #=============================================#
def InitializeFrame(FileNumberArray, DataDirectory, Field, Action):

    
    global time_text
    global elem_text
    global Lines
    global RefLines
    global IC
    global mesh

    # Initialize empty return list
    retlist = []
    
    # Initialize Empty Lines. Add to Return List.
    for i in range(nDirs):
        Lines[i].set_data([],[])
        retlist += [Lines[i]]


    # Initialize time text.  Add to Return List
    time_text.set_text('')
    retlist +=  [ time_text ]

    # Initialize element number text.  Add to Return List
    elem_text.set_text('')
    retlist +=  [ elem_text ]

    
    # If requested, initialize mesh. Add to Return List
    if gvS.PlotMesh:
        mesh.set_data([],[])
        retlist += [ mesh ]
        
        



     # If requested, initialize initial condition lines. Add to Return List
    Data0 = ['None']*nDirs
    DataUnits0 = ['None']*nDirs
    Time0 = ['None']*nDirs
    dX10 = ['None']*nDirs
    X1_C0 = ['None']*nDirs
    if gvS.ShowIC:
        for i in range(nDirs):
            Data0[i], DataUnits0[i],  \
            X1_C0[i], dX10[i], Time0[i] = fetchData_AMReX(0,                 \
                                                          FileNumberArray[i],\
                                                          DataDirectory[i],  \
                                                          Field              )

        Data0, nLines0 = ApplyAction(Data0, Action)
        for i in range(len(Data0)):
            IC[i].set_data( X1_C0[i], Data0[i] )
            retlist += [ IC[i] ]

    # If requested, initialize refinement lines. Add to Return List
    if gvS.ShowRefinement:
        if not gvS.ShowIC:
            for i in range(nDirs):
                Data0[i], DataUnits0[i],  \
                X1_C0[i], dX10[i], Time0[i] = fetchData_AMReX(0,                 \
                                                              FileNumberArray[i],\
                                                              DataDirectory[i],  \
                                                              Field              )
        bottom, top = plt.ylim()
        RefinementLocations = FindRefinementBoundaries( dX10[0] )
        for i in range(len(RefinementLocations)):
            RefLines[i].set_data ( (RefinementLocations[i],  \
                                    RefinementLocations[i] ),\
                                       (top, bottom) )
            retlist +=  [RefLines[i]]


        

    return tuple(retlist)




 #=============================================#
#                                               #
#   UpdateFrame                                 F#
#                                               #
 #=============================================#
def UpdateFrame( t, FileNumberArray, DataDirectory, Field, Action):

    global time_text
    global elem_text
    global Lines
    global RefLines
    global IC
    global mesh


    print('    {:}/{:}'.format( t+1, nSS ) )



    # Draw new value lines.
    retlist = []
    Data = ['None']*nDirs
    DataUnits = ['None']*nDirs
    Time = ['None']*nDirs
    dX1 = ['None']*nDirs
    X1_C = ['None']*nDirs
    for i in range(nDirs):
        Data[i], DataUnits[i],  \
        X1_C[i], dX1[i], Time[i] = fetchData_AMReX(t,                 \
                                                   FileNumberArray[i],\
                                                   DataDirectory[i],  \
                                                   Field              )

    Data, nLines = ApplyAction(Data, Action)

    
    for i in range(nLines):
        Lines[i].set_data( X1_C[i] , Data[i].flatten())
        retlist += [Lines[i]]
        
        
        
    
    # Create new time text.
    time_text.set_text( r'$t={:.3e}\ \left[\mathrm{{{:}}}\right]$' \
                        .format( Time[0], gvU.TimeUnits ) )
    retlist += [ time_text ]

    # Create new element number text.
    elem_text.set_text( r'$Elements: {:}$' \
                        .format( len(X1_C[0] ) ) )
    retlist += [ elem_text ]


    # If requested and amr is true, recreate mesh lines.
    if gvS.PlotMesh and gvS.amr:
        mesh.set_data( X1_C - 0.5 * dX1,            \
                        0.5 * ( vmin + vmax )       \
                          * np.ones( dX1.shape[0] ) )
        retlist += [ mesh ]


    
    # If requested and amr is true, recreate refinement lines.
    if (gvS.ShowRefinement and gvS.amr):
        bottom, top = plt.ylim()
        RefinementLocations = FindRefinementBoundaries( dX1[0] )
        for i in range(len(RefinementLocations)):
            RefLines[i].set_data ( (RefinementLocations[i],  \
                                    RefinementLocations[i] ),\
                                    (top, bottom) )
        retlist += [RefLines[i]]



    return tuple(retlist)



 #=============================================#
#                                               #
#   fetchData_AMReX                             #
#                                               #
 #=============================================#
def fetchData_AMReX(t, FileNumberArray, DataDirectory, Field ):

#    print(FileNumberArray)
#    print(FileNumberArray[t])
    FileDirectory = DataDirectory + str(FileNumberArray[t]) + '/'

    TimeFile = FileDirectory + '{:}.dat'.format( 'Time' )
    X1File   = FileDirectory + '{:}.dat'.format( 'X1' )
    dX1File  = FileDirectory + '{:}.dat'.format( 'dX1' )
    DataFile = FileDirectory + '{:}.dat'.format( Field )

    DataShape, DataUnits, MinVal, MaxVal = ReadHeader( DataFile )

    Time = np.loadtxt( TimeFile )
    X1_C = np.loadtxt( X1File   )
    dX1  = np.loadtxt( dX1File  )
    Data = np.loadtxt( DataFile )

    return Data, DataUnits, X1_C, dX1, Time


 #=============================================#
#                                               #
#   fetchData_Native                             #
#                                               #
 #=============================================#
def fetchData_Native(t, FileNumberArray, DataDirectory, Field ):


    FileDirectory = DataDirectory + str(FileNumberArray[t]) + '/'

    TimeFile = FileDirectory + '{:}.dat'.format( 'Time' )
    X1File   = FileDirectory + '{:}.dat'.format( 'X1' )
    dX1File  = FileDirectory + '{:}.dat'.format( 'dX1' )
    DataFile = FileDirectory + '{:}.dat'.format( Field )

    DataShape, DataUnits, MinVal, MaxVal = ReadHeader( DataFile )

    Time = np.loadtxt( TimeFile )
    X1_C = np.loadtxt( X1File   )
    dX1  = np.loadtxt( dX1File  )
    Data = np.loadtxt( DataFile )

    return Data, DataUnits, X1_C, dX1, Time




 #=============================================#
#                                               #
#   ReadHeader                                  #
#                                               #
 #=============================================#
def ReadHeader( DataFile ):

    f = open( DataFile )

    dum = f.readline()

    s = f.readline(); ind = s.find( ':' )+1
    DataShape = np.array( list( map( np.int64, s[ind:].split() ) ), np.int64 )

    s = f.readline(); ind = s.find( ':' )+1
    DataUnits = s[ind:]

    s = f.readline(); ind = s.find( ':' )+1
    MinVal = np.float64( s[ind:] )

    s = f.readline(); ind = s.find( ':' )+1
    MaxVal = np.float64( s[ind:] )

    f.close()

    return DataShape, DataUnits, MinVal, MaxVal



 #=============================================#
#                                               #
#   ApplyMovieSettings                          #
#                                               #
 #=============================================#
def ApplyMovieSettings( ax, xL, xH, dX10, Field, DataUnits ):

#    global IC
#    global mesh
    
    ax.set_title( r'$\texttt{{{:}}}$'.format( gvS.FigTitle ), fontsize = 15 )
    
    ax.set_xlabel \
    ( r'$x^{{1}}\ \left[\mathrm{{{:}}}\right]$'.format( gvU.X1Units ), fontsize = 15 )


    ax.set_ylabel( Field  + ' ' + r'$\left[\mathrm{{{:}}}\right]$' \
                              .format( DataUnits[2:-2] ) )
    
    ax.legend( prop = {'size':12} )
    ax.grid(which='both')


    ax.set_xlim( xL, xH )
    ax.set_ylim( gvS.vmin, gvS.vmax )

    #
    #   Apply Log Scale
    #
    if gvS.UseLogScale_Y: ax.set_yscale( 'log' )
    
    if gvS.UseLogScale_X:
        xL = max( xL, xL + 0.25 * dX10[0]/9 )
        ax.set_xlim( xL, xH )
        ax.set_xscale( 'log' )


    #
    #   Apply Extra Lines
    #

    ax.set_ylim( gvS.vmin, gvS.vmax )


    return



 #=============================================#
#                                               #
#   ApplyAction                                 #
#                                               #
 #=============================================#
def ApplyAction( Data, Action):

    if len(Data) > 1:

        if ( Action.lower() == 'reldiff' ):
            NewData = [abs(Data[0][:]-Data[1][:])/abs(Data[0][:])]
            nLines = 1
            
        else:
            NewData = Data
            nLines = len(Data)
    else:
        NewData = Data
        nLines = len(Data)

    return NewData, nLines
