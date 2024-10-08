#!/usr/bin/env python3

import numpy as np
from sys import argv
from matplotlib import animation
from functools import partial
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

import GlobalVariables.Settings as gvS
import GlobalVariables.Units    as gvU
import Utilities.BounceFinder   as BF
from Utilities.RefinementBoundaryFinder import FindRefinementBoundaries
from Utilities.Synchronizer import SynchronizeFrameLists
from Utilities.FetchData import fetchData_AMReX, fetchData_Native, ReadHeader


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
    global Frame_List
    global nFrames
    global Data0
    
    nSS = ['None']*nDirs
    for i in range(nDirs):
        if DataDirectory[i][-1] != '/': DataDirectory[i] += '/'
        nSS[i] = len(FileNumberArray[i])
    

    Data0, DataUnits, X1_C0, dX10, Time = fetchData_AMReX(0,                 \
                                                          FileNumberArray[0],\
                                                          DataDirectory[0],  \
                                                          Field              )
                                                          
    

#    Frame_List = SynchronizeFrameLists(nDirs,nSS,'Time')
    Frame_List = SynchronizeFrameLists( nDirs,              \
                                        nSS,                \
                                        FileNumberArray,    \
                                        DataDirectory,      \
                                        'Bounce'            )
        
    nFrames = len(Frame_List)


    if not gvS.UseCustomLimits:
        gvS.vmin = +np.inf
        gvS.vmax = -np.inf
        for i in range( nDirs ):
            for j in range( nSS[i] ):
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
                                    frames = nFrames,                               \
                                    blit = True )



    fps = max( 1, nFrames / gvS.MovieRunTime )
    


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

    ColorList = ['blue','red','green']
    LineList  = ['solid','dashed','dotted']
    LabelList = ['AMR, 0.5 km','Uni, 0.50 km','AMR, 1.00 km']
    
    elem_offset = 0.0
    if gvS.ReferenceBounce:
        elem_offset = 0.1





    Lines = ['None']*nDirs
    time_text = ['None']*nDirs
    elem_text = ['None']*nDirs
    for i in range(nDirs):
        time_text[i] = ax.text( 0.05, 0.20-0.05*i, '', transform = ax.transAxes, fontsize = 13 )
        elem_text[i] = ax.text( 0.45+elem_offset, 0.20-0.05*i, '', transform = ax.transAxes, fontsize = 13 )
        Lines[i], = ax.plot( [],[],                                                         \
                             color     = ColorList[i],                                     \
                             linestyle = LineList[i],                                       \
                             marker    = '',                                                \
                             label     = r'$u_{{{:}}}\left(t\right)$'.format(LabelList[i]),     \
                             zorder    = 10 )
    
    
    if gvS.PlotMesh:
        mesh,     = ax.plot( [],[] )
    
    
    if gvS.ShowIC:
        IC = ['None']*nDirs
        for i in range(nDirs):
            IC[i], = ax.plot( [],[],                                                            \
                              color  = 'k',                                                     \
                              linestyle = LineList[i],                                          \
                              label  = r'$u_{{{:}}}\left(0\right)$'.format(LabelList[i]),     \
                              zorder = 1000,                                                    \
                              alpha  = 0.5    )
    
    
    if gvS.ShowRefinement:
#        RefLines = [['None']*gvS.RefinementLevels]*nDirs
        RefLines = [['None'] * gvS.RefinementLevels for i in range(nDirs)]
            
        for j in range(nDirs):
            for i in range(gvS.RefinementLevels):
                RefLines[j][i], = ax.plot( [],[],                       \
                                           color     = ColorList[j],    \
                                           linestyle = LineList[j],     \
                                           scaley    = False,           \
                                           zorder    = 0,               \
                                           alpha     = 0.7              )


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
    global f_old
    global t_old
    global Data0
    global nDirs

    # Initialize empty return list
    retlist = []
    
    # Initialize Empty Lines. Add to Return List.
    for i in range(nDirs):
        Lines[i].set_data([],[])
        retlist += [Lines[i]]

        # Initialize time text.  Add to Return List
        time_text[i].set_text('')
        retlist +=  [ time_text[i] ]

        # Initialize element number text.  Add to Return List
        elem_text[i].set_text('')
        retlist +=  [ elem_text[i] ]

    
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
    for i in range(nDirs):
        Data0[i], DataUnits0[i],  \
        X1_C0[i], dX10[i], Time0[i] = fetchData_AMReX(0,                 \
                                                      FileNumberArray[i],\
                                                      DataDirectory[i],  \
                                                      Field              )
                                                      

    Data1, nLines0 = ApplyAction(Data0, Field, Action)
    
    for i in range(nLines0):
        Lines[i].set_data( X1_C0[i] , Data1[i].flatten())
        retlist += [Lines[i]]


    if gvS.ShowIC:
        for i in range(nLines0):
            IC[i].set_data( X1_C0[i], Data1[i] )
            retlist += [ IC[i] ]

    # If requested, initialize refinement lines. Add to Return List
    if gvS.ShowRefinement:
        bottom, top = plt.ylim()
        for j in range(nDirs):
            RefinementLocations = FindRefinementBoundaries( dX10[j] )
            for i in range(len(RefinementLocations)):
                RefLines[j][i].set_data ( [RefinementLocations[i],  \
                                           RefinementLocations[i]],\
                                           (top, bottom) )
                retlist +=  [RefLines[j][i]]


    

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
    global f_old
    global t_old
    global nFrames
    global Frame_List
    global Data0
    global nDirs
    

    print('Creating Frame: {:}/{:} \r'.format( t+1, nFrames ), end='' )



    # Draw new value lines.
    retlist = []
    Data = ['None']*nDirs
    DataUnits = ['None']*nDirs
    Time = ['None']*nDirs
    dX1 = ['None']*nDirs
    X1_C = ['None']*nDirs
    
    
    for i in range(nDirs):
        Data[i], DataUnits[i],  \
        X1_C[i], dX1[i], Time[i] = fetchData_AMReX(Frame_List[t][i],  \
                                                   FileNumberArray[i],\
                                                   DataDirectory[i],  \
                                                   Field              )

    Data, nLines = ApplyAction(Data, Field, Action)

    for i in range(nLines):
        Lines[i].set_data( X1_C[i] , Data[i].flatten())
        retlist += [Lines[i]]
        
        # Create new time text.
        if gvS.ReferenceBounce:
            time_text[i].set_text( r'$t_{{{:}}}-t_{{b}}  ={:.5e}\ \left[\mathrm{{{:}}}\right]$' \
                                .format( gvS.DataType[i][0], Time[i]-BF.BounceTimeList[i], gvU.TimeUnits ) )
        else:
            time_text[i].set_text( r'$t_{{{:}}}  ={:.5e}\ \left[\mathrm{{{:}}}\right]$' \
                                .format( gvS.DataType[i][0], Time[i], gvU.TimeUnits ) )
        retlist += [ time_text[i] ]




        # Create new element number text.
        elem_text[i].set_text( r'$N_{{{:}}}={:}$' \
                            .format( gvS.DataType[i][0], len(X1_C[i] ) ) )
        retlist += [ elem_text[i] ]





    # If requested and amr is true, recreate mesh lines.
    if gvS.PlotMesh and gvS.amr:
        mesh.set_data( X1_C - 0.5 * dX1,            \
                        0.5 * ( vmin + vmax )       \
                          * np.ones( dX1.shape[0] ) )
        retlist += [ mesh ]


    
    # If requested and amr is true, recreate refinement lines.
    if (gvS.ShowRefinement and gvS.amr):
        bottom, top = plt.ylim()
        for j in range(nDirs):
            RefinementLocations = FindRefinementBoundaries( dX1[j] )
            for i in range(len(RefinementLocations)):
                RefLines[j][i].set_data ( [RefinementLocations[i],  \
                                           RefinementLocations[i] ],\
                                           (top, bottom) )
                retlist += [RefLines[j][i]]



    return tuple(retlist)







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
    
    ax.legend(  prop = {'size':12},         \
                loc = "upper right"          )
    
    ax.grid(which='both')
#    ax.tight_layout(rect=[0, 0, 0.75, 1])

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
def ApplyAction( Data, Field, Action):

    global Data0

    if len(Data) > 1:

        if ( Action.lower() == 'reldiff' ):
            NewData = [abs(Data[0][:]-Data[1][:])/abs(Data[0][:])]
            nLines = 1
            
        elif (Field == 'PolytropicConstant'):

            NewData = ['None']*nDirs
            for i in range(nDirs):
                NewData[i] = Data[i]/Data0[i]
            nLines = len(NewData)
        else:
            NewData = Data
            nLines = len(Data)
    else:
        NewData = Data
        nLines = len(Data)

    return NewData, nLines
