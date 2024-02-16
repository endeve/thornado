#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )


import GlobalVariables.Settings as gvS
import GlobalVariables.Units    as gvU

from Utilities.GetPlotData  import GetPlotData
from Utilities.Files        import GetFileNumberArray
from Utilities.MakeDataArray import MakeProbelmDataDirectory

import Utilities.DecadeFinder   as DF
import Utilities.TimeFinder     as TF

if __name__ == "__main__":
    #### ========== User Input ==========

    # Specify name of problem
    ProblemName = 'AdiabaticCollapse_XCFC'

    # Specify title of figure
    gvS.FigTitle = ProblemName

    THORNADO_DIR = '/Users/nickroberts/thornado/'
    #THORNADO_DIR = '/home/kkadoogan/Work/Codes/thornado/'

    # Specify directory containing amrex Plotfiles
    gvS.nDirs = 2


    # Specify directory containing amrex Plotfiles
    PlotDirectories = ['None']*gvS.nDirs
    PlotDirectories[0] \
      = THORNADO_DIR + \
          'SandBox/AMReX/Applications/AdiabaticCollapse_XCFC/AdiabaticCollapse_XCFC_AMR_dr0.25km'
          
    PlotDirectories[1] \
      = THORNADO_DIR + \
          'SandBox/AMReX/Applications/AdiabaticCollapse_XCFC/AdiabaticCollapse_XCFC_Uni_dr0.5km'
          
          
    gvS.DataType = ['None']*gvS.nDirs
    gvS.DataType[0] = 'AMReX'
    gvS.DataType[1] = 'AMReX'

    # Specify plot file base name
    PlotBaseName = ProblemName + '.plt'

    # Specify field to plot
    Field = 'PF_D'

    # Specify to plot in log-scale
    gvS.UseLogScale_X  = True
    gvS.UseLogScale_Y  = True
    gvS.UseLogScale_2D = False

    # Specify whether or not to use physical units
    UsePhysicalUnits = True

    # Specify coordinate system (currently supports 'cartesian' and 'spherical')
    CoordinateSystem = 'cartesian'

    # Max level of refinement to plot (-1 plots leaf elements)
    gvS.MaxLevel = -1

    # Only use every <plotEvery> plotfile
    PlotEvery = ['1']*gvS.nDirs
    PlotEvery[0] = 2000
    PlotEvery[1] = 1000
    
    # First and last snapshots and number of snapshots to include in movie
    SSi = -1 # -1 -> SSi = 0
    SSf = -1 # -1 -> plotfileArray.shape[0] - 1
    nSS = -1 # -1 -> plotfileArray.shape[0]
    
    # Write extra info to screen
    gvS.Verbose = True

    # Use custom limts for y-axis (1D) or colorbar (2D)
    gvS.UseCustomLimits = False
    gvS.vmin            = 0.0
    gvS.vmax            = 20.0

    gvS.ShowRefinement = False
    gvS.RefinementLevels = 9

    # Save figure (True) or plot figure (False)
    gvS.SaveFig = False

    # Specify colormap (2D only)
    cmap = 'viridis'

    
    #### ====== End of User Input =======
    
    labelList = ['0.25km AMR','0.5km Uni']
    
    colorlist = ['blue','red','green']
    markerlist = ['solid','dashed','dashdot']
    
    polar = False
    if CoordinateSystem == 'spherical':
        polar = True

    ID      = '{:s}_{:s}'.format( ProblemName, Field )
    FigName = 'fig.{:s}_Many.png'.format( ID )

    
    DataDirectories = ['None']*gvS.nDirs

    DataDirectories[0] = 'DataDirectories/{:s}_AMR025km_AMR'.format( ProblemName )
    DataDirectories[1] = 'DataDirectories/{:s}_AMR05km_Uni'.format( ProblemName )

    # Append "/" to PlotDirectory, if not present
    for i in range(gvS.nDirs):
        if not PlotDirectories[i][-1] == '/': PlotDirectories[i] += '/'


    gvU.SetSpaceTimeUnits(CoordinateSystem, UsePhysicalUnits)


    FileNumberArrays = [['None']]*gvS.nDirs
    FileNumberArraysB = [['None']]*gvS.nDirs
    #   Plot by the density Decade
    Decades, DecadeFrames, DecadeTimes = DF.CreateDecadeData( PlotDirectories,  \
                                                              PlotBaseName,     \
                                                              DataDirectories,  \
                                                              gvS.DataType      )

    FileNumberArrays = DecadeFrames
    

    TimeList = [50, 100, 150, 200, 300, 400] # ms
    TimeFrames = TF.FindTimeFrames( TimeList,           \
                                    PlotDirectories,    \
                                    PlotBaseName,       \
                                    gvS.DataType        )
                                    
    FileNumberArraysB = TimeFrames

#   PlotEvery between SSi and SSf
#    for i in range(gvS.nDirs):
#        FileNumberArrays[i] = GetFileNumberArray( PlotDirectories[i],     \
#                                                  PlotBaseName,       \
#                                                  SSi, SSf,           \
#                                                  PlotEvery[i]           )


    print(FileNumberArrays)
    print(FileNumberArraysB)
    
    exit()
    fig, ax = plt.subplots( 1, 1 )


    nDirs = len(FileNumberArrays)
    
    for i in range(nDirs):
        nFiles = len(FileNumberArrays[i])

        for f in range(nFiles):

            Data, DataUnit, Time,   \
            X1_C, X2_C, X3_C,       \
            dX1, dX2, dX3 = GetPlotData( PlotDirectories[i],     \
                                         PlotBaseName,      \
                                         Field,             \
                                         argv = ['a',str(FileNumberArrays[i][f])],   \
                                         Verbose = gvS.Verbose,     \
                                         DataType = gvS.DataType[i]  )



            nX1 = X1_C.shape[0]
            nX2 = X2_C.shape[0]
            nX3 = X3_C.shape[0]


            nDims = 1
            if nX2 > 1: nDims += 1
            if nX3 > 1: nDims += 1

            if nDims == 1:

                ax.plot( X1_C, Data,        \
                         scaley = False,    \
                         color  = colorlist[i],    \
                         zorder = 0,        \
                         alpha  = 0.4,      \
                         label  = labelList[i], \
                         linestyle = markerlist[i]       )
                
            else:

                print('Multi-D not implemented yet.')
                
             
    xL = X1_C[0 ] - 0.5 * dX1[0 ]
    xH = X1_C[-1] + 0.5 * dX1[-1]
         
    if gvS.UseLogScale_X:
        ax.set_xscale( 'log' )
        xL = max( xL, 0.0 + 0.25 * (X1_C[1]-X1_C[0]) )
        
    if gvS.UseLogScale_Y: ax.set_yscale( 'log' )

    if gvS.UseCustomLimits: ax.set_ylim( vmin, vmax )

    if gvS.ShowRefinement:
        bottom, top = plt.ylim()
        ax.plot( (gvS.RefinementLocations[:], gvS.RefinementLocations[:]), \
                 (top, bottom),     \
                 scaley = False,    \
                 color  = 'red',    \
                 zorder = 0,        \
                 alpha  = 0.4       )


    ax.set_xlim( xL, xH )
    ax.set_xlabel \
      ( r'$x^{{1}}\ \left[\mathrm{{{:}}}\right]$'.format( gvU.X1Units ), \
        fontsize = 15 )
        
    ax.set_ylabel( Field + ' ' + '$'+DataUnit+'$' )
    ax.grid()


    ax.legend(  prop = {'size':12},         \
            loc = "upper right"          )





    if gvS.SaveFig:

        plt.savefig( FigName, dpi = 300 )

    else:

        plt.show()

    plt.close()

    import os
    os.system( 'rm -rf __pycache__ ' )

